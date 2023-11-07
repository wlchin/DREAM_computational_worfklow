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
smoldata_dotplot_select = smoladata[smoladata.obs["predicted.celltype.l2"].isin(celltypes_selected)]
fraction = smoldata_dotplot_select.obs.shape[0] * 0.20
smoller = sc.pp.filter_genes(smoldata_dotplot_select, min_counts=None, min_cells=fraction, max_counts=None, max_cells=None, inplace=False, copy=True)

def plot_signatures_on_common_celltypes(gene_signatureR, gene_signatureNR, savefile): 
    blood_responder = pd.read_csv(gene_signatureR, header = None)[0]
    blood_nonresponder = pd.read_csv(gene_signatureNR, header = None)[0]
    markers_express_blood_responder = adata.var[smoller[0]].gene[adata.var[smoller[0]].gene.isin(blood_responder)]
    markers_express_blood_nonresponder = adata.var[smoller[0]].gene[adata.var[smoller[0]].gene.isin(blood_nonresponder)]
    #markers_express_tumour_responder = adata.var[smoller[0]].gene[adata.var[smoller[0]].gene.isin(tumour_responder)]

    #res_only = ranks_df[ranks_df["features"].isin(markers_express_blood_responder)]
    #non_res_only = ranks_df[ranks_df["features"].isin(markers_express_blood_nonresponder)]

    marker_dict = {"top in R" : markers_express_blood_responder, 
                  "top in NR" : markers_express_blood_nonresponder}

    sc.pl.dotplot(smoldata_dotplot_select, marker_dict, 
                  groupby='predicted.celltype.l2', 
                  dendrogram=False,
                 standard_scale="var", swap_axes=True, save=savefile)


plot_signatures_on_common_celltypes("data/blood_responder_sampled_25genes.txt", "data/blood_nonresponder_sampled_25genes.txt", "blood_markers_factor_loadings_25_genes.png")
plot_signatures_on_common_celltypes("data/blood_responder_sampled_50genes.txt", "data/blood_nonresponder_sampled_50genes.txt", "blood_markers_factor_loadings_50_genes.png")
plot_signatures_on_common_celltypes("data/blood_responder_sampled_100genes.txt", "data/blood_nonresponder_sampled_100genes.txt", "blood_markers_factor_loadings_100_genes.png")

# plot CD8 T cells

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
res = pd.read_csv("data/blood_responder_sampled_50genes.txt", header = None)[0]

markers_express_blood_responder = CD8TEM.var["gene"][CD8TEM.var.index.isin(res)].to_list()

#markers_express_blood_nonresponder = pd.read_csv("results/blood_markers_for_dotplot_NR.txt", header = None)[0]

sc.pl.dotplot(CD8TEM, markers_express_blood_responder, 
            groupby='Mapped_Category', 
            dendrogram=False,
            standard_scale="var", swap_axes=True, save="CD8TEM_blood_markers.png")