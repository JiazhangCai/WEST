from SpaceFlow import SpaceFlow
import scanpy as sc
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata
import copy
import random
import torch
import SpaGCN as spg
import warnings
import datatable as dt
warnings.filterwarnings("ignore")
from numpy import array
from pyWNN import pyWNN
from sklearn.manifold import MDS
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import OneHotEncoder

class west(object):

    def __init__(self, 
                X, # expression matrix
                S): # spatial information
        super(west, self).__init__()
        S.names = ["x", "y"]
        self.X = X 
        self.S = S  
        adata = anndata.AnnData(pd.DataFrame(X.to_numpy().astype(float)))
        adata.obs["loc_x"] =  S["x"].to_numpy()
        adata.obs["loc_y"] =  S["y"].to_numpy()
        self.adata = adata

    def preprocess(self, 
                    min_cells = 3, # number of minimum expressed spots for each gene 
                    n_top_genes = 3000, # number of remained genes with highest number of counts
                    n_comps = 100, # number of remained components in PCA
                    screening = True): # whether implements screening to genes
        adata = self.adata
        adata.var_names_make_unique()
        spg.prefilter_genes(adata, min_cells=min_cells)    # avoiding all genes are zeros
        spg.prefilter_specialgenes(adata)    # only keep the genes whose name is started with "ERCC" or "MT-"
        sc.pp.normalize_per_cell(adata)    # Normalization
        sc.pp.log1p(adata)    # Take log
        if screening:
            sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, subset=True)    # Gene screening
        sc.pp.pca(adata, n_comps=n_comps)    # PCA
        self.adata = adata 
        self.min_cells = min_cells
        self.n_top_genes = n_top_genes
        self.n_comps = n_comps

    def train(self, 
                n_clusters = 10, # number of clusters for SpaGCN
                resolution = 0.3, # resolution for SpaceFlow
                num_pcs = [100,100], # embedded dimensions for SpaGCN and SpaceFlow
                n_neighbors = 100, # number of neighbors for the weighted ensemble method
                p = 0.5, # spatial factor for SpaGCN
                regularization = 0.01, # spatial factor for SpaceFlow
                seed = 1234): # random seed
        adata = self.adata
        X = self.X
        S = self.S

        # set seed
        random.seed(seed)
        torch.manual_seed(seed)
        np.random.seed(seed)

        # SpaGCN
        adata_SpaGCN = adata
        adj = spg.pairwise_distance(S.to_numpy().astype(np.float32))     # calculate distance matrix
        l = spg.search_l(p, adj, start=0.001, end=1000, tol=0.01, max_run=100)     # calculate l based on p
        res = spg.search_res(adata_SpaGCN, adj, l, n_clusters, start=0.3, step=0.1, tol=5e-3, lr=0.05, max_epochs=200, r_seed=seed, t_seed=seed, n_seed=seed)    # calculate resolution based on l

        # train
        clf = spg.SpaGCN()
        clf.set_l(l)
        clf.train(adata_SpaGCN, adj, num_pcs=num_pcs[0], init_spa=True, init="louvain", res=res, tol=5e-3, lr=0.05, max_epochs=500)

        y_pred, prob=clf.predict()

        adata.obs["SpaGCN"]= y_pred
        adata.obs["SpaGCN"]=adata.obs["SpaGCN"].astype("category")

        # get new embedding
        z,q = clf.model.predict(clf.embed,clf.adj_exp)
        adata.obsm["SpaGCN"] = copy.deepcopy(z.detach().numpy())

        # SpaceFlow
        sf = SpaceFlow(expr_data=pd.DataFrame(X.to_numpy()), spatial_locs=S.to_numpy())
        sf.preprocessing_data(n_top_genes=self.n_top_genes, flavor='seurat', n_comps=self.n_comps, n_neighbors=n_neighbors)

        # train
        sf.train(spatial_regularization_strength=regularization, 
                z_dim=num_pcs[1], 
                lr=0.001, 
                epochs=1000, 
                max_patience=100, 
                min_stop=100, 
                random_seed=seed, 
                gpu=0)

        sf.segmentation(domain_label_save_filepath="./domains1.tsv", 
                n_neighbors=n_neighbors, 
                resolution=resolution)

        est_label = np.array(sf.domains).astype(int)
        adata.obs["SpaceFlow"] = est_label
        adata.obs["SpaceFlow"] = adata.obs["SpaceFlow"].astype("category")

        # get new embedding
        adata.obsm["SpaceFlow"] = sf.embedding

        # WNN
        WNNobj = pyWNN(adata, reps=["SpaceFlow", "SpaGCN"], npcs=num_pcs, n_neighbors=n_neighbors, seed=seed)
        adata = WNNobj.compute_wnn(adata)
        
        self.adata = adata

    def cluster(self, 
                resolution=1): # resolution for clustering 
        adata = self.adata
        sc.tl.leiden(adata, key_added="WNN", obsp="WNN", resolution=resolution)
        self.adata = adata
    
    def embedding(self, 
                dim=50): # dimension of embedding
        similarity_matrix = np.array(self.adata.obsp["WNN"].toarray())
        distance_matrix = similarity_matrix + similarity_matrix.T
        mds = MDS(n_components=dim, dissimilarity='precomputed')
        embedding = mds.fit_transform(distance_matrix)
        
        return embedding
