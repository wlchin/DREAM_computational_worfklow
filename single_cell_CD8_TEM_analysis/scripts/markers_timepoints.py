import scanpy as sc
import pandas as pd
import gseapy as gs
import numpy as np
from scipy.stats import zscore

print("am reading")

CD8 = sc.read("data/CD8TEM.h5ad")
CD8.var.index = CD8.var["gene"]
clusters = pd.read_csv("data/CD8_seurat_clusters.csv")

clusters.index = clusters["Unnamed: 0"]
CD8 = CD8[clusters["Unnamed: 0"]]
CD8.obs["seurat_clusters"] = clusters["cluster_seurat"]

CD8_clean = CD8[~((CD8.obs["seurat_clusters"] == 3) | (CD8.obs["seurat_clusters"] == 6))]

sc.pp.normalize_per_cell(CD8_clean, counts_per_cell_after=1e4)
sc.pp.log1p(CD8_clean)

CD8_clean.obs["seurat_clusters"] = CD8_clean.obs["seurat_clusters"].astype("category")

#sc.tl.rank_genes_groups(CD8_clean, 'seurat_clusters', 
def get_cluster_markers(CD8_clean, clusterlab, savefile):
    CD8_TEM_2 = CD8_clean[CD8_clean.obs["seurat_clusters"] == clusterlab]
    CD8_TEM_2.obs["PFS_6M_cat"] = CD8_TEM_2.obs["PFS_6M_Timepoint"].astype("category")
    res = pd.DataFrame(columns=CD8_TEM_2.var_names, index=CD8_TEM_2.obs['PFS_6M_cat'].cat.categories)      

    for clust in CD8_TEM_2.obs.PFS_6M_cat.cat.categories: 
        res.loc[clust] = CD8_TEM_2[CD8_TEM_2.obs['PFS_6M_Timepoint'].isin([clust]),:].X.mean(0)
    df = res[["CD81", "HLA-A", "HLA-DRB1", "CCL5", "CD6", "EOMES", "HLA-E", "AKT1", "HLA-DQA1", "IRF1", "RHOH", "RUNX3", "CD8B", "IFNG", "CD5"]]
    """
    df = res[['CD6', 'CCL4', 'HLA-A', 'CD81', 
                   'CD300A', 'RPS6', 'HLA-E', 'FCRL3', 'PRNP', 'PKN1',
                  'XCL1', 'SMAD7',
             'BTN3A2', 'PLCB1', 'RPS3', 'TNF', 'CD160', 'FGR', 'HSPB1', 'CD244', 'FLOT1', 'SASH3','PLCG2']] # terminally exhausted
    """
    df = df.astype("float32")
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    df_z = df[numeric_cols].apply(zscore)
    df_z.to_csv(savefile)

#GO:0042110 T cell activation

get_cluster_markers(CD8_clean, 0, "results/res_cluster0.csv")
get_cluster_markers(CD8_clean, 1, "results/res_cluster1.csv")
get_cluster_markers(CD8_clean, 2, "results/res_cluster2.csv")