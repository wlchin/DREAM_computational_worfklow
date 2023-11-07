import pandas as pd
import matplotlib.pyplot as plt

import scanpy as sc
import pandas as pd
import anndata as ad

sc.set_figure_params(scanpy=True, dpi_save=400)

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

adata = ad.read("data/dream.h5ad")
adata.var.index = adata.var.gene
smoladata = adata[adata.obs["percent.mt"] < 10]

smoladata = adata[adata.obs["percent.mt"] < 10]
smoldata_dotplot_select = smoladata[smoladata.obs["predicted.celltype.l2"].isin(celltypes_selected)]
fraction = smoldata_dotplot_select.obs.shape[0] * 0.20
tumour_responder = pd.read_csv("data/tumour_responder_sampled_50genes.txt", header = None)[0]

smoller = sc.pp.filter_genes(smoldata_dotplot_select, min_counts=None, min_cells=fraction, max_counts=None, max_cells=None, inplace=False, copy=True)
markers_express_tumour_responder = adata.var[smoller[0]].gene[adata.var[smoller[0]].gene.isin(tumour_responder)]
#markers_express_tumour_responder = adata.var[smoller[0]].gene[adata.var[smoller[0]].gene.isin(tumour_responder)]

#res_only = ranks_df[ranks_df["features"].isin(markers_express_blood_responder)]
#non_res_only = ranks_df[ranks_df["features"].isin(markers_express_blood_nonresponder)]

sc.pl.dotplot(smoldata_dotplot_select, markers_express_tumour_responder, 
              groupby='predicted.celltype.l2', 
              dendrogram=False,
             standard_scale="var", swap_axes=True, save="tumour_genes_50.png")
