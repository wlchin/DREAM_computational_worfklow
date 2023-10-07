import scanpy as sc
import pandas as pd
import anndata as ad
import pickle

### this is the CD8 TEM in seurat 

adata = ad.read("data/dream.h5ad")
adata.var.index = adata.var.gene
smoladata = adata[adata.obs["percent.mt"] < 10]

CD8_seurat_clusterlabels = pd.read_csv("data/CD8_seurat_clusters.csv")
indices_of_CD8_T_cells = CD8_seurat_clusterlabels[CD8_seurat_clusterlabels["cluster_seurat"].isin([0,1,2])]["Unnamed: 0"]

#indices_of_CD8_T_cells = CD8_seurat_clusterlabels["Unnamed: 0"]

CD8TEM = smoladata[smoladata.obs.index.isin(indices_of_CD8_T_cells)]
smolCD8TEM = CD8_seurat_clusterlabels[CD8_seurat_clusterlabels["cluster_seurat"].isin([0,1,2])]
smolCD8TEM.index = smolCD8TEM["Unnamed: 0"]
CD8TEM.obs["cluster_seurat"] = smolCD8TEM["cluster_seurat"].astype("string")
label_mapping = {'0': 'CD8_TEM_1', '1': 'CD8_TEM_2', '2': 'CD8_TEM_3'}

# Use the map() function to create a new column with mapped labels
CD8TEM.obs['Mapped_Category'] = CD8TEM.obs["cluster_seurat"].map(label_mapping).astype("category")

#fraction = smoldata_dotplot_select.obs.shape[0] * 

#tumour_responder = pd.read_csv("tumour_responder.txt", header = None)[0]
#smoller = sc.pp.filter_genes(smoldata_dotplot_select, min_counts=None, min_cells=fraction, max_counts=None, max_cells=None, inplace=False, copy=False)
#markers_express_tumour_responder = adata.var[smoller[0]].gene[adata.var[smoller[0]].gene.isin(tumour_responder)]

markers_express_tumour_responder = pd.read_csv("results/tumour_markers_for_dotplot.txt", header = None)[0]
markers_express_blood_responder = pd.read_csv("results/blood_markers_for_dotplot.txt", header = None)[0]


sc.pl.dotplot(CD8TEM, markers_express_blood_responder, 
              groupby='Mapped_Category', 
              dendrogram=False,
             standard_scale="var", swap_axes=True, save="CD8TEM_blood_markers.png")

sc.pl.dotplot(CD8TEM, markers_express_tumour_responder, 
              groupby='Mapped_Category', 
              dendrogram=False,
             standard_scale="var", swap_axes=True, save="CD8TEM_tumour_markers.png")