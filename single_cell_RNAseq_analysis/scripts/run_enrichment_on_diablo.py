import scanpy as sc
import pandas as pd
import anndata as ad
import decoupler as dc
from scipy.stats import zscore

adata = ad.read("data/dream.h5ad")

adata.var.index = adata.var.gene
smoladata = adata[adata.obs["percent.mt"] < 10]
CD8TEM = smoladata[smoladata.obs["predicted.celltype.l2"] == "CD8 TEM"]

def perform_enrichment(signature_file, savefile_common_celltypes, savefile_CD8TEM):
    R_sig = pd.read_csv(signature_file, header = None)
    R_sig = R_sig.drop_duplicates()
    R_sig["geneset"] = "DIABLO"
    R_sig.columns = ["genesymbol", "geneset"]
    dc.run_ora(mat=smoladata, net=R_sig, source='geneset', target='genesymbol', verbose=True, use_raw=False, min_n=20)
    acts = dc.get_acts(smoladata, obsm_key='ora_estimate')
    GSEA = acts.obsm["ora_estimate"]
    GSEA["celltype"] = adata.obs["predicted.celltype.l2"]
    GSEA["Timepoint"] = adata.obs["Timepoint"]
    df_expression = GSEA.groupby(["celltype", "Timepoint"]).mean().unstack().drop_duplicates().dropna()
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
    df_expression.loc[celltypes_selected].to_csv(savefile_common_celltypes)

    dc.run_ora(mat=CD8TEM, net=R_sig, source='geneset', target='genesymbol', verbose=True, use_raw=False, min_n=20)
    acts = dc.get_acts(CD8TEM, obsm_key='ora_estimate')
    GSEA = acts.obsm["ora_estimate"]
    GSEA["celltype"] = adata.obs["predicted.celltype.l2"]
    GSEA["Timepoint"] = adata.obs["Timepoint"]
    CD8_seurat_clusterlabels = pd.read_csv("data/CD8_seurat_clusters.csv")
    CD8_seurat_clusterlabels.columns = ["source", "seurat_cluster"]
    CD8_seurat_clusterlabels.index = CD8_seurat_clusterlabels.source
    totCD8 = pd.concat([CD8_seurat_clusterlabels, GSEA], axis = 1, join="inner")
    totCD8[["seurat_cluster", "DIABLO", "Timepoint"]]
    df_expression_CD8 = totCD8.groupby(["seurat_cluster", "Timepoint"]).mean().unstack().drop_duplicates().dropna()
    df_expression_CD8.to_csv(savefile_CD8TEM)
    
perform_enrichment("data/blood_responder_sampled_25genes.txt", "results/GSEA_DIABLO_25_genes.csv", "results/GSEA_DIABLO_CD8_25_genes.csv")
perform_enrichment("data/blood_responder_sampled_50genes.txt", "results/GSEA_DIABLO_50_genes.csv", "results/GSEA_DIABLO_CD8_50_genes.csv")
perform_enrichment("data/blood_responder_sampled_100genes.txt", "results/GSEA_DIABLO_100_genes.csv", "results/GSEA_DIABLO_CD8_100_genes.csv")
