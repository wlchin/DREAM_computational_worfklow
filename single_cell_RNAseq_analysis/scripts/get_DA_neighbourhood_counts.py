import numpy as np
import pandas as pd
import scanpy as sc
from tqdm import tqdm as tqdm

adata_neighbourhoods = sc.read(snakemake.input[0])

def get_single_neighbourhood(adata_neighbourhoods, ind_of_interest):
    ind = adata_neighbourhoods.obsm["nhoods"][:,ind_of_interest].toarray().ravel()
    vals = adata_neighbourhoods.obs.iloc[np.where(ind)].groupby("Patient").size()
    neighbourhood_name = "index-" + str(ind_of_interest)
    x = (vals/sum(ind)).to_frame(neighbourhood_name)
    return x

def get_multiple_neighbourhoods(adata_neighbourhoods, index_list):
    index_res = []
    for i in tqdm(index_list):
        a = get_single_neighbourhood(adata_neighbourhoods, i)
        index_res.append(a)
    total = pd.concat(index_res, axis = 1).fillna(0)
    return total

def get_single_neighbourhood_counts(adata_neighbourhoods, ind_of_interest):
    ind = adata_neighbourhoods.obsm["nhoods"][:,ind_of_interest].toarray().ravel()
    vals = adata_neighbourhoods.obs.iloc[np.where(ind)].groupby("Patient").size()
    neighbourhood_name = "index-" + str(ind_of_interest)
    x = vals.to_frame(neighbourhood_name)
    return x

def get_multiple_neighbourhoods_counts(adata_neighbourhoods, index_list):
    index_res = []
    for i in tqdm(index_list):
        a = get_single_neighbourhood_counts(adata_neighbourhoods, i)
        index_res.append(a)
    total = pd.concat(index_res, axis = 1).fillna(0)
    return total

## timepoint 0 

FDR_cutoff = 0.1

all_index = pd.read_csv(snakemake.input[1])
indicer = all_index[all_index["FDR"] < FDR_cutoff]["Unnamed: 0"]
test = get_multiple_neighbourhoods_counts(adata_neighbourhoods, indicer.to_list())

test.columns = all_index[all_index["FDR"] < FDR_cutoff]["index_cell"]
test.index = ["Pt_" + str(i) for i in test.index.to_list()]
test.transpose().to_csv(snakemake.output[0])

## timepoint 1

all_index = pd.read_csv(snakemake.input[2])
indicer = all_index[all_index["FDR"] < FDR_cutoff]["Unnamed: 0"]
test = get_multiple_neighbourhoods_counts(adata_neighbourhoods, indicer.to_list())

test.columns = all_index[all_index["FDR"] < FDR_cutoff]["index_cell"]
test.index = ["Pt_" + str(i) for i in test.index.to_list()]
test.transpose().to_csv(snakemake.output[1])

## timepoint 2

all_index = pd.read_csv(snakemake.input[3])
indicer = all_index[all_index["FDR"] < FDR_cutoff]["Unnamed: 0"]
test = get_multiple_neighbourhoods_counts(adata_neighbourhoods, indicer.to_list())

test.columns = all_index[all_index["FDR"] < FDR_cutoff]["index_cell"]
test.index = ["Pt_" + str(i) for i in test.index.to_list()]
test.transpose().to_csv(snakemake.output[2])
