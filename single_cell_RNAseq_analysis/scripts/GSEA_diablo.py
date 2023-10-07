import scanpy as sc
import pandas as pd
import anndata as ad
import decoupler as dc
from scipy.stats import zscore

adata = ad.read("data/dream.h5ad")

adata.var.index = adata.var.gene

genesig = pd.read_csv("data/test_blood_tp12.txt", header = None)
genesig2 = pd.read_csv("data/test_blood2_tp12.txt", header = None)

genesig["geneset"] = "DIABLO_tp1"
genesig2["geneset"] = "DIABLO_tp1"

combined_df = pd.concat([genesig, genesig2], axis = 0)

stuff = combined_df.drop_duplicates()

stuff.columns = ["genesymbol", "geneset"]

smoladata = adata[adata.obs["percent.mt"] < 10]

dc.run_ora(mat=smoladata, net=stuff, source='geneset', target='genesymbol', verbose=True, use_raw=False, min_n=20)

acts = dc.get_acts(smoladata, obsm_key='ora_estimate')

GSEA = acts.obsm["ora_estimate"]

GSEA["celltype"] = adata.obs["predicted.celltype.l2"]
GSEA["Timepoint"] = adata.obs["Timepoint"]

df_expression = GSEA.groupby(["celltype", "Timepoint"]).mean().unstack().drop_duplicates().dropna()

zscore_df = df_expression.apply(zscore)


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

df_expression.loc[celltypes_selected].to_csv("results/GSEA_DIABLO.csv")

# get CD8 T cells and seurat clusters

CD8_seurat_clusterlabels = pd.read_csv("data/CD8_seurat_clusters.csv")
CD8_seurat_clusterlabels.columns = ["source", "seurat_cluster"]
CD8_seurat_clusterlabels.index = CD8_seurat_clusterlabels.source
#dc.run_ora(mat=CD8TEM, net=stuff, source='geneset', target='genesymbol', verbose=True, use_raw=False, min_n=20)

totCD8 = pd.concat([CD8_seurat_clusterlabels, GSEA], axis = 1, join="inner")
totCD8[["seurat_cluster", "DIABLO_tp1", "Timepoint"]]
df_expression_CD8 = totCD8.groupby(["seurat_cluster", "Timepoint"]).mean().unstack().drop_duplicates().dropna()
df_expression_CD8.to_csv("results/GSEA_DIABLO_CD8.csv")