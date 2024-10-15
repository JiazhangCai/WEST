import scanpy as sc
import numpy as np
import copy
from sklearn import preprocessing
from scipy.sparse import csr_matrix, lil_matrix, diags
import time
from sklearn.manifold import MDS


def get_nearestneighbor(knn, neighbor=1):
    '''For each row of knn, returns the column with the lowest value
    I.e. the nearest neighbor'''
    indices = knn.indices
    indptr = knn.indptr
    data = knn.data
    nn_idx = []
    for i in range(knn.shape[0]):
        cols = indices[indptr[i]:indptr[i+1]]
        rowvals = data[indptr[i]:indptr[i+1]]
        idx = np.argsort(rowvals)
        nn_idx.append(cols[idx[neighbor-1]])
    return(np.array(nn_idx))


def compute_bw(knn_adj, embedding, n_neighbors=20):
    intersect = knn_adj.dot(knn_adj.T)
    indices = intersect.indices
    indptr = intersect.indptr
    data = intersect.data
    data = data / ((n_neighbors*2) - data)
    bandwidth = []
    for i in range(intersect.shape[0]):
        cols = indices[indptr[i]:indptr[i+1]]
        rowvals = data[indptr[i]:indptr[i+1]]
        idx = np.argsort(rowvals)
        valssort = rowvals[idx]
        numinset = len(cols)
        if numinset<n_neighbors:
            sys.exit('Fewer than 20 cells with Jacard sim > 0')
        else:
            curval = valssort[n_neighbors]
            for num in range(n_neighbors, numinset):
                if valssort[num]!=curval:
                    break
                else:
                    num+=1
            minjacinset = cols[idx][:num]
            if num <n_neighbors:
                print('shouldnt end up here')
                sys.exit(-1)
            else:
                euc_dist = ((embedding[minjacinset,:]-embedding[i,:])**2).sum(axis=1)**.5
                euc_dist_sorted = np.sort(euc_dist)[::-1]
                bandwidth.append( np.mean(euc_dist_sorted[:n_neighbors]) )
    return(np.array(bandwidth))



def compute_affinity(dist_to_predict, dist_to_nn, bw):
    affinity = dist_to_predict-dist_to_nn
    affinity[affinity<0]=0
    affinity = affinity * -1
    affinity = np.exp(affinity / (bw-dist_to_nn))
    return(affinity)


def dist_from_adj(adjacency, embed, nndist):
    result = []
    indices = adjacency.indices
    indptr = adjacency.indptr
    ncells = adjacency.shape[0]
    
    tic = time.perf_counter()
    for k in range(len(embed)):
        dist = lil_matrix(adjacency.shape)
        for i in range(ncells):
            for j in range(indptr[i], indptr[i+1]):
                col = indices[j]
                a = (((embed[k][i,:] - embed[k][col,:])**2).sum()**.5) - nndist[k][i]
                if a == 0: dist[i,col] = np.nan
                else: dist[i,col] = a
        result.append(dist) 
                
    toc = time.perf_counter()
    print('time: %.2f seconds' % (toc-tic))

    return(result)


def select_topK(dist,  n_neighbors=20):
    indices = dist.indices
    indptr = dist.indptr
    data = dist.data
    nrows = dist.shape[0]

    final_data = []
    final_col_ind = []

    tic = time.perf_counter()
    for i in range(nrows):
        cols = indices[indptr[i]:indptr[i+1]]
        rowvals = data[indptr[i]:indptr[i+1]]
        idx = np.argsort(rowvals)
        final_data.append(rowvals[idx[(-1*n_neighbors):]])
        final_col_ind.append(cols[idx[(-1*n_neighbors):]])
            
    final_data = np.concatenate(final_data)
    final_col_ind = np.concatenate(final_col_ind)
    final_row_ind = np.tile(np.arange(nrows), (n_neighbors, 1)).reshape(-1, order='F')
                
    result = csr_matrix((final_data, (final_row_ind, final_col_ind)), shape=(nrows, dist.shape[1]))

    return(result)


class WEST():
    
    def __init__(self, adata, reps, n_neighbors=20, npcs=None, seed=1234, distances=None):
        
        """
        :param adata: an anndata.AnnData type object
        :param reps: the embedding to integrate, e.g., ['SpaGCN', 'SpaceFlow']
        :param n_neighbors: number of neighbors, defalut: 20
        :param npcs: number of PCs to use, e.g., [20, 30]. If None, then use the original embedding, defalut: None
        :param seed: random seed, default: 1234
        :param distances: distance matrix. If None, then calculate the distance matrix, defalut: None
        """
        
        self.seed = seed
        np.random.seed(seed)
        
        self.adata = adata.copy()
        self.reps = reps
        if npcs is None:
            self.npcs = [self.adata.obsm[key].shape[1] for key in reps]
        else:
            self.npcs = npcs
            
        for (i,r) in enumerate(reps):
            self.adata.obsm[self.reps[i]] = preprocessing.normalize(adata.obsm[r][:,0:self.npcs[i]])

        self.n_neighbors = n_neighbors
        
        #### Calculate distance if not provided
        self.distances = []
        if distances is None:
            print('Computing KNN distance matrices using default Scanpy implementation')
            for i in range(len(self.npcs)):
                sc.pp.neighbors(self.adata, n_neighbors=n_neighbors, n_pcs=self.npcs[i], use_rep=self.reps[i], metric='euclidean', key_added=self.reps[i])    
                self.distances.append(self.reps[i]+'_distances')
            
            
        else:
            self.distances = distances      
        
        #### Convert to sparse matrix 
        for d in self.distances:
            if type(self.adata.obsp[d]) is not csr_matrix:
                self.adata.obsp[d] = csr_matrix(self.adata.obsp[d])
            
        self.NNdist = [] # distance to the nearest neighbor
        self.NNidx = [] # index of the nearest neighbor
        self.NNadjacency = [] # adjacency matrix
        self.BWs = [] # bandwidth

        for (i,r) in enumerate(self.reps):
            nn = get_nearestneighbor(self.adata.obsp[self.distances[i]])
            dist_to_nn = ((self.adata.obsm[r]-self.adata.obsm[r][nn, :])**2).sum(axis=1)**.5
            nn_adj = (self.adata.obsp[self.distances[i]]>0).astype(int)
            nn_adj_wdiag = nn_adj.copy()
            nn_adj_wdiag.setdiag(1)
            bw = compute_bw(nn_adj_wdiag, self.adata.obsm[r], n_neighbors=self.n_neighbors)
            self.NNidx.append(nn)
            self.NNdist.append(dist_to_nn)
            self.NNadjacency.append(nn_adj)
            self.BWs.append(bw)

        self.weights = []
        self.WEST = None
    
    def compute_weights(self):
        cmap = {i: [j for j in range(len(self.reps)) if j != i] for i in range(len(self.reps))}
        affinity_ratios = []
        for (i,r) in enumerate(self.reps):
            within_predict = self.NNadjacency[i].dot(self.adata.obsm[r]) / (self.n_neighbors-1)
            within_predict_dist = ((self.adata.obsm[r] - within_predict)**2).sum(axis=1)**.5
            within_affinity = compute_affinity(within_predict_dist, self.NNdist[i], self.BWs[i])
            
            affinity_ratio_sub = []
            for j in range(len(cmap[i])):
                cross_predict = self.NNadjacency[cmap[i][j]].dot(self.adata.obsm[r]) / (self.n_neighbors-1)
                cross_predict_dist = ((self.adata.obsm[r] - cross_predict)**2).sum(axis=1)**.5
                cross_affinity = compute_affinity(cross_predict_dist, self.NNdist[i], self.BWs[i])
                affinity_ratio_sub.append(np.exp(within_affinity-(cross_affinity+0.001)))
                
            affinity_ratios.append(sum(affinity_ratio_sub))         
        
        self.weights = [x/sum(affinity_ratios) for x in affinity_ratios]
   
    def compute_west(self):
        print('Computing modality weights')
        self.compute_weights()
        union_adj_mat = (sum(self.adata.obsp[self.distances[i]] for i in range(len(self.reps))) > 0).astype(int)
        
        embed = [self.adata.obsm[self.reps[i]] for i in range(len(self.reps))]
        nndist = self.NNdist
        full_dists = dist_from_adj(union_adj_mat, embed, nndist)
        
        weighted_dist = csr_matrix(union_adj_mat.shape)
        for (i,dist) in enumerate(full_dists):
            dist = diags(-1 / (self.BWs[i] - self.NNdist[i]), format='csr').dot(dist)
            dist.data = np.exp(dist.data)
            ind = np.isnan(dist.data)
            dist.data[ind] = 1
            dist = diags(self.weights[i]).dot(dist)
            weighted_dist += dist

        print('Selecting top K neighbors')
        self.WEST = select_topK(weighted_dist,  n_neighbors=self.n_neighbors)
        WESTdist = self.WEST.copy()
        x = (1-WESTdist.data) / 2
        x[x<0]=0
        x[x>1]=1
        WESTdist.data = np.sqrt(x)
        self.WESTdist = WESTdist

        self.adata.obsp['WEST'] = self.WEST
        self.adata.obsp['WEST_distance'] = self.WESTdist
        for i in range(len(self.reps)):
            self.adata.obsm[self.reps[i]] = self.adata.obsm[self.reps[i]]
        self.adata.uns['WEST'] = {'connectivities_key': 'WEST',
                             'distances_key': 'WEST_distance',
                             'params': {'n_neighbors': self.n_neighbors,
                             'method': 'WEST',
                             'random_state': self.seed,
                             'metric': 'euclidean',
                             'use_rep': self.reps[0],
                             'n_pcs': self.npcs[0]}}
    
    def get_embedding(self, dim=50): 
        similarity_matrix = np.array(self.adata.obsp["WEST"].toarray())
        distance_matrix = (similarity_matrix + similarity_matrix.T) / 2
        mds = MDS(n_components=dim, dissimilarity='precomputed')
        new_embedding = mds.fit_transform(distance_matrix)
        
        return(new_embedding)
    