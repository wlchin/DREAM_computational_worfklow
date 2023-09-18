import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd

blood2 = sc.read("results/pca.h5ad")

col_scheme = pd.read_csv("results/colscheme.csv")
colscheme = dict(zip(col_scheme["celltype"],col_scheme["colour"]))

with plt.rc_context({"figure.figsize": (6,6), "figure.dpi": (400)}):
    sc.pl.umap(blood2, color = ["predicted.celltype.l2"], legend_loc='on data', legend_fontsize = 7, size = 1, palette = colscheme, legend_fontweight = "normal")
    plt.savefig('results/by_celltype.png')
    
blood2.obs["Patient_cat"] = blood2.obs["Patient"].astype("category")
blood2.obs["PFS_cat"] = blood2.obs["PFS_6M"].astype("category")
blood2.obs["Batch_cat"] = blood2.obs["Batch"].astype("category")
blood2.obs["Block_cat"] = blood2.obs["Block"].astype("category")
blood2.obs["Timepoint_cat"] = blood2.obs["Timepoint"].astype("category")

with plt.rc_context({"figure.figsize": (6,6), "figure.dpi": (400)}):
    sc.pl.umap(blood2, color = ["PFS_cat"])
    plt.savefig('results/by_pfs.png')
    
with plt.rc_context({"figure.figsize": (6,6), "figure.dpi": (400)}):
    sc.pl.umap(blood2, color = ["Timepoint_cat"])
    plt.savefig('results/by_timepoint.png')