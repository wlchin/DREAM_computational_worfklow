import scanpy as sc
import pandas as pd
import gseapy as gs
import numpy as np
from scipy.stats import zscore

CD8 = sc.read(snakemake.input[0])
CD8.var.index = CD8.var["gene"]
clusters = pd.read_csv("data/CD8_seurat_clusters.csv")

clusters.index = clusters["Unnamed: 0"]
CD8 = CD8[clusters["Unnamed: 0"]]
CD8.obs["seurat_clusters"] = clusters["cluster_seurat"].astype("category")

CD8_clean = CD8[~((CD8.obs["seurat_clusters"] == 3) | (CD8.obs["seurat_clusters"] == 6))]

sc.pp.normalize_per_cell(CD8_clean, counts_per_cell_after=1e4)
sc.pp.log1p(CD8_clean)

CD8_clean.obs["seurat_clusters"] = CD8_clean.obs["seurat_clusters"].astype("category")

CD8_clean.obs.to_csv("data/metadata_CD8.csv")

#sc.pp.scale(CD8_clean)

adata = CD8_clean

res = pd.DataFrame(columns=adata.var_names, index=adata.obs['seurat_clusters'].cat.categories)                                                                                                 

for clust in adata.obs.seurat_clusters.cat.categories: 
    res.loc[clust] = adata[adata.obs['seurat_clusters'].isin([clust]),:].X.mean(0)

#output = res[["PDCD1", "TOX", "HAVCR2", "LAG3", "TIGIT", "TCF7", "IL7R", "LEF1", "CCR7", "PRF1", "GZMA"]]
#output = res[["PDCD1", "TOX", "HAVCR2", "LAG3", "TIGIT", "TCF7", "IL7R", "LEF1", "CCR7"]]

#df = res[["PDCD1", "TOX", "HAVCR2", "LAG3", "TIGIT", 
#              "TCF7", "IL7R", "LEF1", "CCR7", 
#              "PRF1", "GZMA", "GZMB", "GZMK", "GNLY", "FAS",
#             "SELL", "BATF", "ID3", "BCL6", "NR4A3",
#             "CX3CR1", "CXCR6", "CCL3", "ID2", "ZEB2", "TOX2", "KLRG1", "EOMES",
#             "KIR3DL1", "KIR3DL3", "TYROBP", "KLRC3", "KLRC2", "NCR1", "NCR3", "IKZF2",
#             "CD69", "PTPRC", "TBR1"]] # terminally exhausted

df = res[["PDCD1", "HAVCR2", "LAG3", "TIGIT", "TNFRSF9", "TNFRSF18", "ENTPD1", # immune checkpoint
              "PRF1", "GZMA", "GZMB", "GZMK", "GNLY", "FAS", "ITGAE", #effector
              "CCL3", "CX3CR1", "CXCL13", "CXCR6", "KLRG1", "ID2", "ZEB2", "EOMES", "PRDM1", "TOX", # exhaustion
              "TCF7", "MYB", "SLAMF6", "IL7R", "CCR7", "LEF1", "SELL", "BATF", # stem like
             "KLRC3", "KLRC1", "KLRC2", "KIR3DL1", "KIR3DL3", "NCR1"]] # KLR

# CCL3, CX3CR1, CXCL13,CXCR6, KLRG1, ID2, ZEB2, EOMES, PRDM1, TOX. 

#KeyError: "['41BB', 'GITR', 'CD39', 'CD103', 'CD62L'] not in index"
#GITR TNFRSF18
#41BB TNFRSF9
#CD39 ENTPD1
#CD103 ITGAE
# CD62L is SELL

#IKZF2 and NK cell–associated genes (e.g., TYROBP, KLRC2, KLRC3, NCR1, and NCR3

df.index = ["CD8_TEM_" + str(i) for i in range(1,8)]

df = df.astype("float32")

numeric_cols = df.select_dtypes(include=[np.number]).columns
df_z = df[numeric_cols].apply(zscore)

df_z.to_csv("results/zscaled_proprotions.csv")

print(df_z.shape)