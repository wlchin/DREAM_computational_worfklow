import pandas as pd
import anndata as ad
import numpy as np

adata = ad.read("data/dream.h5ad")

adata.var.index = adata.var.gene

genesig = pd.read_csv("data/test_blood_tp12.txt", header = None)
genesig2 = pd.read_csv("data/test_blood2_tp12.txt", header = None)

genesig["geneset"] = "DIABLO_tp1"
genesig2["geneset"] = "DIABLO_tp1"

combined_df = pd.concat([genesig, genesig2], axis = 0)
combined_df_no_duplicates = combined_df.drop_duplicates()
combined_df_no_duplicates.columns = ["genesymbol", "geneset"]

smoladata = adata[adata.obs["percent.mt"] < 10]

def grouped_obs_mean(adata, group_key, layer=None, gene_symbols=None):
    if layer is not None:
        getX = lambda x: x.layers[layer]
    else:
        getX = lambda x: x.X
    if gene_symbols is not None:
        new_idx = adata.var[idx]
    else:
        new_idx = adata.var_names

    grouped = adata.obs.groupby(group_key)
    out = pd.DataFrame(
        np.zeros((adata.shape[1], len(grouped)), dtype=np.float64),
        columns=list(grouped.groups.keys()),
        index=adata.var_names
    )

    for group, idx in grouped.indices.items():
        X = getX(adata[idx])
        out[group] = np.ravel(X.mean(axis=0, dtype=np.float64))
    return out

df_average = grouped_obs_mean(smoladata, "predicted.celltype.l2")


celltypes_selected = ["CD14 Mono",
               "CD4 TCM",
               "NK",
               "CD8 TEM",
               "CD4 Naive",
               "CD16 Mono",
               "B naive",
               "CD8 Naive",
               "CD4 TEM",
               "Treg"]

ind = df_average.index.isin(combined_df_no_duplicates["genesymbol"])

df_average[ind][celltypes_selected].to_csv("results/average_expression_individual_genes_DIABLO.csv")