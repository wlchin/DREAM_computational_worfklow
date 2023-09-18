import pandas as pd
import scanpy as sc

enrichment = pd.read_csv("results/enrichment_all_celltypes_top_pathways.csv")
sc_data = sc.read("results/condensed.h5ad")

contrast0 = pd.read_csv("results/contrast_0.csv")
contrast1 = pd.read_csv("results/contrast_1.csv")
contrast2 = pd.read_csv("results/contrast_2.csv")

print("timepoint 0")

sc_data.obs.index = contrast0["index_cell"]
contrast0.index = contrast0["index_cell"]
combined_df = pd.concat([sc_data.obs, contrast0], axis = 1, join = "inner")

enrichment.index = combined_df.index
enrichment["celltype"] = combined_df["nhood_annotation"]
enrichment["logFC"] = combined_df["logFC"]
enrichment["FDR"] = combined_df["FDR"]

for i in enrichment["celltype"].unique():
    x = enrichment[enrichment["celltype"] == i].sort_values("logFC")
    filename = "results/" + i.replace(" ", "_") + "_timepoint0.csv"
    x.drop(columns = "Unnamed: 0", inplace = True)
    x.to_csv(filename)
    
print("timepoint 1")

sc_data.obs.index = contrast1["index_cell"]
contrast1.index = contrast1["index_cell"]
combined_df = pd.concat([sc_data.obs, contrast1], axis = 1, join = "inner")

enrichment.index = combined_df.index
enrichment["celltype"] = combined_df["nhood_annotation"]
enrichment["logFC"] = combined_df["logFC"]
enrichment["FDR"] = combined_df["FDR"]

for i in enrichment["celltype"].unique():
    x = enrichment[enrichment["celltype"] == i].sort_values("logFC")
    filename = "results/" + i.replace(" ", "_") + "_timepoint1.csv"
    x.drop(columns = "Unnamed: 0", inplace = True)
    x.to_csv(filename)
    
### timepoint 2

sc_data.obs.index = contrast2["index_cell"]
contrast2.index = contrast2["index_cell"]
combined_df = pd.concat([sc_data.obs, contrast2], axis = 1, join = "inner")

enrichment.index = combined_df.index
enrichment["celltype"] = combined_df["nhood_annotation"]
enrichment["logFC"] = combined_df["logFC"]
enrichment["FDR"] = combined_df["FDR"]

for i in enrichment["celltype"].unique():
    x = enrichment[enrichment["celltype"] == i].sort_values("logFC")
    filename = "results/" + i.replace(" ", "_") + "_timepoint2.csv"
    x.drop(columns = "Unnamed: 0", inplace = True)
    x.to_csv(filename)