import scanpy as sc
import pandas as pd
import anndata as ad
import decoupler as dc


adata = ad.read("data/dream.h5ad")
adata.var.index = adata.var.gene
smoladata = adata[adata.obs["percent.mt"] < 10]

CD8TEM = smoladata[smoladata.obs["predicted.celltype.l2"] == "CD8 TEM"]

print(CD8TEM.shape)

combined_df = pd.read_csv("data/blood_responder.txt", header = None)
combined_df["geneset"] = "DIABLO_blood"
stuff = combined_df.drop_duplicates()

print(stuff.head())

stuff.columns = ["genesymbol", "geneset"]

dc.run_ora(mat=CD8TEM, net=stuff, source='geneset', target='genesymbol', verbose=True, use_raw=False, min_n=20)
acts = dc.get_acts(CD8TEM, obsm_key='ora_estimate')

GSEA = acts.obsm["ora_estimate"]
GSEA["celltype"] = adata.obs["predicted.celltype.l2"]
GSEA["Timepoint"] = adata.obs["Timepoint"]

CD8_seurat_clusterlabels = pd.read_csv("data/CD8_seurat_clusters.csv")
CD8_seurat_clusterlabels.columns = ["source", "seurat_cluster"]
CD8_seurat_clusterlabels.index = CD8_seurat_clusterlabels.source

totCD8 = pd.concat([CD8_seurat_clusterlabels, GSEA], axis = 1, join="inner")
totCD8[["seurat_cluster", "DIABLO_blood", "Timepoint"]]
df_expression_CD8 = totCD8.groupby(["seurat_cluster", "Timepoint"]).mean().unstack().drop_duplicates().dropna()
df_expression_CD8.to_csv("results/GSEA_DIABLO_CD8.csv")