import pandas as pd
import scanpy as sc

enrichment = pd.read_csv(snakemake.input[0])
sc_data = sc.read(snakemake.input[1])
contrast0 = pd.read_csv(snakemake.input[2])
contrast1 = pd.read_csv(snakemake.input[3])
contrast2 = pd.read_csv(snakemake.input[4])

### timepoint 0 

sc_data.obs["index_cell"] = sc_data.obs.index
combined_df = pd.concat([sc_data.obs, contrast0], axis = 1, join = "inner")

enrichment.index = combined_df.index
enrichment["celltype"] = combined_df["predicted.celltype.l2"]
enrichment["logFC"] = combined_df["logFC"]
enrichment["FDR"] = combined_df["FDR"]

for i in enrichment["celltype"].unique():
    x = enrichment[enrichment["celltype"] == i].sort_values("logFC")
    filename = "results/" + i.replace(" ", "_") + "_timepoint0.csv"
    x.drop(columns = "Unnamed: 0", inplace = True)
    x.to_csv(filename)
    
### timepoint 1

sc_data.obs["index_cell"] = sc_data.obs.index
combined_df = pd.concat([sc_data.obs, contrast1], axis = 1, join = "inner")

enrichment.index = combined_df.index
enrichment["celltype"] = combined_df["predicted.celltype.l2"]
enrichment["logFC"] = combined_df["logFC"]
enrichment["FDR"] = combined_df["FDR"]

for i in enrichment["celltype"].unique():
    x = enrichment[enrichment["celltype"] == i].sort_values("logFC")
    filename = "results/" + i.replace(" ", "_") + "_timepoint1.csv"
    x.drop(columns = "Unnamed: 0", inplace = True)
    x.to_csv(filename)
    
### timepoint 2

sc_data.obs["index_cell"] = sc_data.obs.index
combined_df = pd.concat([sc_data.obs, contrast2], axis = 1, join = "inner")

enrichment.index = combined_df.index
enrichment["celltype"] = combined_df["predicted.celltype.l2"]
enrichment["logFC"] = combined_df["logFC"]
enrichment["FDR"] = combined_df["FDR"]

for i in enrichment["celltype"].unique():
    x = enrichment[enrichment["celltype"] == i].sort_values("logFC")
    filename = "results/" + i.replace(" ", "_") + "_timepoint2.csv"
    x.drop(columns = "Unnamed: 0", inplace = True)
    x.to_csv(filename)