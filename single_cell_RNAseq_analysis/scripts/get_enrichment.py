import scanpy as sc
import pandas as pd
import decoupler as dc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import numpy as np

x = sc.read('data/dream.h5ad')
x = x[~(x.obs["predicted.celltype.l2"] == "Doublet")]
x = x[x.obs["predicted.celltype.l2"] == "CD8 TEM"]
x.var.index = x.var.gene

msigdb = pd.read_pickle("data/msigdb.p")
dc.run_ora(mat=x, net=msigdb, source='geneset', target='genesymbol', verbose=True, use_raw=False)
acts = dc.get_acts(x, obsm_key='ora_estimate')

GSEA = acts.obsm["ora_estimate"]

GSEA["celltype"] = x.obs["predicted.celltype.l2"]
GSEA["PFS_6M"] = x.obs["PFS_6M"]
GSEA["Timepoint"] = x.obs["Timepoint"]

CD8_TEM = GSEA[["celltype", 
                "PFS_6M", 
                "Timepoint", 
                "HALLMARK_INTERFERON_GAMMA_RESPONSE",]]

CD8_TEM[(CD8_TEM["celltype"] == "CD8 TEM")].to_csv("results/CD8TEM_IFNgamma.csv")