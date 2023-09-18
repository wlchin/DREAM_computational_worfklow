from scipy.sparse import csr_matrix
import anndata as ad
import pandas as pd
import scanpy as sc
import milopy
import milopy.core as milo

blood2 = sc.read(snakemake.input[0])

blood2.obsm["X_ref.umap"] = blood2.obsm["X_umap"]
blood2.obs["X_coord"] = blood2.obsm["X_ref.umap"][:,1]
blood2.obs["Y_coord"] = blood2.obsm["X_ref.umap"][:,0]

blood2.obs["milo_sample"] = blood2.obs["Patient"].astype(str) + "_" + blood2.obs["Timepoint"].astype(str)
blood2.obs["contrast_test"] =  "level" + "_" + blood2.obs["PFS_6M_Timepoint"].astype("str")

milo.make_nhoods(blood2)
milo.count_nhoods(blood2, sample_col="milo_sample")

blood2 = blood2.raw.to_adata()
blood2.write(snakemake.output[1])

average_count_mat = csr_matrix(blood2.X.T)
nhoods_X = average_count_mat.dot(blood2.obsm["nhoods"])
denom = blood2.obsm["nhoods"].toarray().sum(0)
tots = nhoods_X/denom

adata_neighbourhoods = ad.AnnData(tots.T)
adata_neighbourhoods.var.index = blood2.var.index
adata_neighbourhoods.var["gene"] = blood2.var.gene
adata_neighbourhoods.var_names = adata_neighbourhoods.var["gene"]

adata_neighbourhoods.write(snakemake.output[0])

pfs6m = blood2.obs[["PFS_6M", "Patient"]].drop_duplicates()
pfs6m["Pt_index"] = "Pt_" + pfs6m["Patient"].astype("string")
pfs6m.index = pfs6m["Pt_index"]
pfs6m.to_pickle("data/index.p")