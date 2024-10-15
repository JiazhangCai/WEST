import pandas as pd
import anndata
import matplotlib.pyplot as plt
import scanpy as sc
import umap


def WEST_scatter(adata, reps='WEST', title=None, cmap=None, s=5, width=5, height=5, switch_coordinate=False, invert_x=False, invert_y=False, save_path=None):
    adata_plt = adata.copy()
    adata_plt.obs[reps] = pd.Categorical(adata_plt.obs[reps])

    #### Define color map if color map is not given
    if cmap is None:
        plot_color = ["#F56867", "#FEB915", "#C798EE", "#59BE86", "#7495D3", "#D1D1D1",
                    "#6D1A9C", "#15821E", "#3A84E6", "#997273", "#787878", "#DB4C6C",
                    "#9E7A7A", "#554236", "#AF5F3C", "#93796C", "#F9BD3F", "#DAB370",
                    "#877F6C", "#268785"]

    #### Get unique categories
    unique_categories = adata_plt.obs[reps].unique()

    #### Create a dictionary mapping categories to colors
    color_dict = {cat: plot_color[i % len(plot_color)] for i, cat in enumerate(unique_categories)}

    #### Set figure size
    plt.figure(figsize=(width, height))

    #### Plot each category separately to assign colors
    if switch_coordinate:
        for category in unique_categories:
            subset = adata_plt.obs[adata_plt.obs[reps] == category]
            plt.scatter(subset['loc_y'], subset['loc_x'], c=color_dict[category], label=category, s=s, alpha=1)
    else:
        for category in unique_categories:
            subset = adata_plt.obs[adata_plt.obs[reps] == category]
            plt.scatter(subset['loc_x'], subset['loc_y'], c=color_dict[category], label=category, s=s, alpha=1)
        
    #### Invert axis
    if invert_x:
        plt.gca().invert_xaxis()
    if invert_y:
        plt.gca().invert_yaxis()

    #### Add a legend
    legend = plt.legend(title='Clusters', bbox_to_anchor=(1.05, 1), loc='upper left')
    frame = legend.get_frame()
    frame.set_edgecolor('black')  # Set the edge color of the legend box to black
    frame.set_linewidth(1.0)  # Set the linewidth of the legend box edge

    #### Add labels and title
    plt.xlabel('x')
    plt.ylabel('y')
    if title is not None:
        plt.title(title)

    #### Save the plot
    if save_path is not None:
        plt.savefig(save_path, bbox_inches='tight')
    
    #### Show the plot
    plt.show()



def WEST_umap(adata, reps='WEST', title=None, cmap=None, s=5, width=5, height=5, save_path=None):
    adata_plt = adata.copy()
    adata_plt.obs[reps] = pd.Categorical(adata_plt.obs[reps])

    #### Define color map if color map is not given
    if cmap is None:
        plot_color = ["#F56867", "#FEB915", "#C798EE", "#59BE86", "#7495D3", "#D1D1D1",
                    "#6D1A9C", "#15821E", "#3A84E6", "#997273", "#787878", "#DB4C6C",
                    "#9E7A7A", "#554236", "#AF5F3C", "#93796C", "#F9BD3F", "#DAB370",
                    "#877F6C", "#268785"]

    #### Get unique categories
    unique_categories = adata_plt.obs['true cluster'].unique()

    #### Create a dictionary mapping categories to colors
    color_dict = {cat: plot_color[i % len(plot_color)] for i, cat in enumerate(unique_categories)}

    #### Compute UMAP coordinates
    if reps=='WEST':
        sc.tl.umap(adata_plt, neighbors_key='WEST')
    else:
        sc.pp.neighbors(adata_plt, use_rep=reps)
        sc.tl.umap(adata_plt)

    #### Set figure size
    plt.figure(figsize=(width, height))

    #### Plot each category separately to assign colors 
    for category in unique_categories:
        subset = adata_plt.obs[adata_plt.obs['true cluster'] == category].index
        subset_indices = [adata_plt.obs.index.get_loc(i) for i in subset]
        plt.scatter(adata_plt.obsm['X_umap'][subset_indices, 0], adata_plt.obsm['X_umap'][subset_indices, 1], c=color_dict[category], label=category, s=s, alpha=1)

    #### Add a legend
    legend = plt.legend(title='Clusters', bbox_to_anchor=(1.05, 1), loc='upper left')
    frame = legend.get_frame()
    frame.set_edgecolor('black')  # Set the edge color of the legend box to black
    frame.set_linewidth(1.0)  # Set the linewidth of the legend box edge

    #### Add labels and title
    plt.xlabel('UMAP1')
    plt.ylabel('UMAP2')
    if title is not None:
        plt.title(title)

    #### Save the plot
    if save_path is not None:
        plt.savefig(save_path, bbox_inches='tight')
    
    #### Show the plot
    plt.show()
