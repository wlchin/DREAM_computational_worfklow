import pandas as pd
import matplotlib.pyplot as plt

import scanpy as sc
import pandas as pd
import anndata as ad

sc.set_figure_params(scanpy=True, dpi_save=400, figsize = (3, 8))

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

ranks_df = pd.read_csv("data/top_rank_loadings_df_blood.csv")
R_genes = pd.read_csv("data/blood_responder_sampled_25genes.txt", header = None)
NR_genes = pd.read_csv("data/blood_nonresponder_sampled_25genes.txt", header = None)

smoldata_dotplot_select = smoladata[smoladata.obs["predicted.celltype.l2"].isin(celltypes_selected)]
fraction = smoldata_dotplot_select.obs.shape[0] * 0.15
blood_responder = R_genes[0]
blood_nonresponder = NR_genes[0]

smoller = sc.pp.filter_genes(smoldata_dotplot_select, min_counts=None, min_cells=fraction, max_counts=None, max_cells=None, inplace=False, copy=True)
markers_express_blood_responder = adata.var[smoller[0]].gene[adata.var[smoller[0]].gene.isin(blood_responder)]
markers_express_blood_nonresponder = adata.var[smoller[0]].gene[adata.var[smoller[0]].gene.isin(blood_nonresponder)]
#markers_express_tumour_responder = adata.var[smoller[0]].gene[adata.var[smoller[0]].gene.isin(tumour_responder)]


res_only = ranks_df[ranks_df["features"].isin(markers_express_blood_responder)]
non_res_only = ranks_df[ranks_df["features"].isin(markers_express_blood_nonresponder)]

res_only = res_only.groupby('features')['weights'].max().reset_index()
non_res_only = non_res_only.groupby('features')['weights'].max().reset_index()

combined_df = pd.concat([res_only, non_res_only], axis = 0)
combined_df["weights"] = combined_df["weights"] * -1
sorted_df = combined_df.sort_values("weights", ascending=False)

finalmarkersNR = sorted_df[sorted_df["weights"] < 0].features
finalmarkersR = sorted_df[sorted_df["weights"] > 0].features

marker_dict = {"R - top genes" : finalmarkersR, 
              "NR - top genes" : finalmarkersNR}

sc.pl.dotplot(smoldata_dotplot_select, marker_dict, 
              groupby='predicted.celltype.l2', 
              dendrogram=False,
             standard_scale="var", swap_axes=True, save="blood_markers_factor_loadings.png")

fig, ax = plt.subplots(figsize = (3,8))

sorted_df = sorted_df.sort_values(by='weights', ascending=False)

categories = sorted_df["features"].to_list()
values = sorted_df["weights"].to_list()
# Create the barplot

data = list(zip(categories, values))
sorted_data = sorted(data, key=lambda x: x[1])
sorted_categories, sorted_values = zip(*sorted_data)
bars = ax.barh(sorted_categories, sorted_values, color=['red' if v > 0 else 'blue' for v in sorted_values])

#bars = ax.barh(sorted_categories, sorted_values)


ax.tick_params(left=False, bottom=False)  # Remove tickmarks
ax.xaxis.set_ticks_position('none')      # Remove x-axis ticks
ax.yaxis.set_ticks_position('none')      # Remove y-axis ticks
ax.grid(False)                          # Remove gridlines

# Add labels and title
ax.set_xlabel('Importance scores')
ax.set_ylabel('Genes')
ax.set_title('Feature importance')

for spine in ax.spines.values():
    spine.set_visible(False)

plt.savefig('results/barplot_feature_importance.png', bbox_inches='tight', dpi = 600)