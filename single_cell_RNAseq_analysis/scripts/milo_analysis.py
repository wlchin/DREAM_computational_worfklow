import scanpy as sc
import numpy as np
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

# types of contrasts
milo.DA_nhoods(blood2, design="~contrast_test", model_contrasts = "contrast_testlevel_4 - contrast_testlevel_3 - contrast_testlevel_1 + contrast_testlevel_0")
milopy.utils.annotate_nhoods(blood2, anno_col='predicted.celltype.l2')
blood2.uns["nhood_adata"].obs.to_csv("results/contrast_0to1.csv")

milo.DA_nhoods(blood2, design="~contrast_test", model_contrasts = "contrast_testlevel_5 - contrast_testlevel_4 - contrast_testlevel_2 + contrast_testlevel_1")
milopy.utils.annotate_nhoods(blood2, anno_col='predicted.celltype.l2')
blood2.uns["nhood_adata"].obs.to_csv("results/contrast_1to2.csv")

milo.DA_nhoods(blood2, design="~contrast_test", model_contrasts = "contrast_testlevel_3 - contrast_testlevel_0") #0
milopy.utils.annotate_nhoods(blood2, anno_col='predicted.celltype.l2')
blood2.uns["nhood_adata"].obs.to_csv("results/contrast_0.csv")

milo.DA_nhoods(blood2, design="~contrast_test", model_contrasts = "contrast_testlevel_4 - contrast_testlevel_1") #1
milopy.utils.annotate_nhoods(blood2, anno_col='predicted.celltype.l2')
blood2.uns["nhood_adata"].obs.to_csv("results/contrast_1.csv")

milo.DA_nhoods(blood2, design="~contrast_test", model_contrasts = "contrast_testlevel_5 - contrast_testlevel_2") #2
milopy.utils.annotate_nhoods(blood2, anno_col='predicted.celltype.l2')
blood2.uns["nhood_adata"].obs.to_csv("results/contrast_2.csv")

