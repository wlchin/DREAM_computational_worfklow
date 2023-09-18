import scanpy as sc

blood = sc.read(snakemake.input[0]) # this has index correct

CD8TEM = blood[blood.obs["predicted.celltype.l2"] == "CD8 TEM"]

CD8TEM.write("data/CD8_TEM.h5ad")
CD8TEM.obs.to_csv("data/metadata_CD8TEM.csv")