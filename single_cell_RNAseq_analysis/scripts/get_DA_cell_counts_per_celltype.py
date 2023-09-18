import pandas as pd
import scanpy as sc
import seaborn as sns

lfc = pd.read_csv(snakemake.input[0])

DA_neighbourhoods = lfc[lfc["FDR"] < 0.1]
binned_neighbhourhoods_on_LFC = pd.cut(DA_neighbourhoods["logFC"], [-5,-4, -3, -2, -1, 0, 1, 2, 3, 4, 5])

patients = pd.read_csv(snakemake.input[1])
LFC_and_FDR_dataframe = DA_neighbourhoods[["Unnamed: 0", "logFC", "FDR"]]
LFC_and_FDR_dataframe.columns = ["index_cell", "logFC", "FDR"]
ptslfc = patients.merge(LFC_and_FDR_dataframe)

binned_neighbhourhoods_on_LFC2 = pd.cut(ptslfc["logFC"], [-5,-4, -3, -2, -1, 0, 1, 2, 3, 4, 5])
patientmat = ptslfc.drop(['logFC', "FDR"], axis = 1).groupby(binned_neighbhourhoods_on_LFC2).sum()
pheno = pd.read_pickle("data/index.p").reindex(patientmat.T.index)

colour_code = colours = {0:"r", 1:"b"}
row_colours = pheno["PFS_6M"].map(colour_code)

fig = sns.clustermap(patientmat.T, row_colors=row_colours, col_cluster=False)
fig.savefig(snakemake.output[0])

csv_for_vis_heatmap = patientmat.T
csv_for_vis_heatmap["pheno"] = pheno["PFS_6M"]
csv_for_vis_heatmap.to_csv(snakemake.output[1])