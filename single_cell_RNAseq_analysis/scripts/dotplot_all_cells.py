import scanpy as sc
import pandas as pd
import anndata as ad

adata = ad.read("data/dream.h5ad")
adata.var.index = adata.var.gene
smoladata = adata[adata.obs["percent.mt"] < 10]

celltypes_selected = ["CD14_Mono",
               "CD4_TCM",
               "NK",
               "CD8_TEM",
               "CD4_Naive",
               "CD16_Mono",
               "B_naive",
               "CD8_Naive",
               "CD4_TEM",
               "Treg"]

smoldata_dotplot_select = smoladata[smoladata.obs["predicted.celltype.l2"].isin(celltypes_selected)]

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

fraction = smoldata_dotplot_select.obs.shape[0] * 0.5
smoldata_dotplot_select = smoladata[smoladata.obs["predicted.celltype.l2"].isin(celltypes_selected)]
blood_responder = pd.read_csv("data/blood_responder.txt", header = None)[0]
blood_nonresponder = pd.read_csv("data/blood_nonresponder.txt", header = None)[0]
tumour_responder = pd.read_csv("data/tumour_responder.txt", header = None)[0]

smoller = sc.pp.filter_genes(smoldata_dotplot_select, min_counts=None, min_cells=fraction, max_counts=None, max_cells=None, inplace=False, copy=False)
markers_express_blood_responder = adata.var[smoller[0]].gene[adata.var[smoller[0]].gene.isin(blood_responder)]
markers_express_blood_nonresponder = adata.var[smoller[0]].gene[adata.var[smoller[0]].gene.isin(blood_nonresponder)]
markers_express_tumour_responder = adata.var[smoller[0]].gene[adata.var[smoller[0]].gene.isin(tumour_responder)]

marker_dict = {"top genes in R" : markers_express_blood_responder,
              "top genes in NR" : markers_express_blood_nonresponder}


sc.pl.dotplot(smoldata_dotplot_select, marker_dict, 
              groupby='predicted.celltype.l2', 
              dendrogram=False,
             standard_scale="var", swap_axes=True, save="blood_markers.png")

fraction = smoldata_dotplot_select.obs.shape[0] * 0.1
tumour_responder = pd.read_csv("data/tumour_responder.txt", header = None)[0]
smoller = sc.pp.filter_genes(smoldata_dotplot_select, min_counts=None, min_cells=fraction, max_counts=None, max_cells=None, inplace=False, copy=False)
markers_express_tumour_responder = adata.var[smoller[0]].gene[adata.var[smoller[0]].gene.isin(tumour_responder)]

sc.pl.dotplot(smoldata_dotplot_select, markers_express_tumour_responder, 
              groupby='predicted.celltype.l2', 
              dendrogram=False,
             standard_scale="var", swap_axes=True, save="tumour_markers.png")

with open("results/blood_markers_for_dotplot.txt", 'w') as file:
    for item in markers_express_blood_responder:
        file.write(item + '\n')

with open("results/tumour_markers_for_dotplot.txt", 'w') as file:
    for item in markers_express_tumour_responder:
        file.write(item + '\n')
