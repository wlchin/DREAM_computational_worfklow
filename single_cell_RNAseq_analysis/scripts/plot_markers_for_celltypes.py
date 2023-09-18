import scanpy as sc
import numpy as np
import pandas as pd


sc.set_figure_params(dpi_save=600, format = 'png', figsize = (10,15))

adata = sc.read("data/dream.h5ad") 
adata.var_names = adata.var["gene"]

adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

sc.pp.filter_cells(adata, min_genes=200)
blood = adata[adata.obs.pct_counts_mt < 10, :]
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)


MonoDC = {'CD16+ monocyte': ["FCGR3A"],
 'CD14+ monocyte': ['CD14', "LYZ", "S100A8"],
 'Plasmacytoid dendritic cell(pDC)': ['CD80', 'CD86', 'CD40', 'CD83'],
 'Conventional dendritic cell(cDC1)': ['CLEC9A'],
 'Conventional dendritic cell(cDC2)': ["CD1C", "CD33"]}

cell_type = ['CD14 Mono',
 'CD16 Mono',
 'cDC2',
 'cDC1',
 'pDC', 
 "ASDC"]

mono = adata[adata.obs["predicted.celltype.l2"].isin(cell_type)]
sc.pl.dotplot(mono, MonoDC, 'predicted.celltype.l2', dendrogram=False, save="Mono")

cell_type = [
 'CD8 TEM',
 'CD8 Proliferating',
 'CD8 Naive',
 'CD8 TCM']

CD8 = {'Naive CD8+ T cell': ['TCF7', 'CCR7'],
 'Central memory CD8+ T cell': ['CD27', 'SELL', "SOCS3"],
 'Exhausted CD8+ T cell': ['TIGIT', 'HAVCR2', 'LAG3', 'PDCD1'],
 'Terminal effector CD8+ T cell': ['NKG7', 'GZMA', 'PRF1', 'GZMB'],
 'Activated CD8+ T cell': ['CD8A', 'CCL5', 'CD3D', 'HLA-DRA']}

cd8 = adata[adata.obs["predicted.celltype.l2"].isin(cell_type)]
sc.pl.dotplot(cd8, CD8, 'predicted.celltype.l2', dendrogram=False, save="CD8")

CD4 = {'Naive CD4 T cell': ['TCF7', 'MAL', 'CCR7', 'MYC', 'TSHZ2'],
 'Memory CD4 T cell': ["CD69", "ANXA1", "S1PR1", "CCR7"],
 'Activated CD4+ T cell': ["IL2RA", "CD38", "CXCR3", "TFRC"],
 'CD4+ cytotoxic T cell': ['FGFBP2']}

cell_type = ['CD4 TCM',
 'CD4 Naive',
 'CD4 Proliferating',
 'CD4 CTL',
 'CD4 TEM']

cd4 = adata[adata.obs["predicted.celltype.l2"].isin(cell_type)]
sc.pl.dotplot(cd4, CD4, 'predicted.celltype.l2', dendrogram=False, save="CD4")

cell_type = [
 'B intermediate',
 'B memory',
 'B naive',
 'Plasmablast']

BCELL = {'Naive B cell': ['CD38', 'CD19', 'CD24', 'CD27', 'IGHD'],
         "Activated B cell" :["CD27", 'CD19'],
        'Memory B cell': ["AICDA", "CD86", "FCER1G"]}

bcell = adata[adata.obs["predicted.celltype.l2"].isin(cell_type)]
sc.pl.dotplot(bcell, BCELL, 'predicted.celltype.l2', dendrogram=False, save="Bcell")


cell_type =  ['NK_CD56bright',
 'NK Proliferating',
 'NK']

NK = {'Natural killer cell': ['FCGR3A', "NCR1"],
     "NK-CD56-Bright":["NCAM1", "XCL1"],
     "NK-CD56-Dim":["IL21R", "GZMB", "KIR3DL1"],
     "NKT":["CD3E", "NKG7", "CD8A"]}

nk = adata[adata.obs["predicted.celltype.l2"].isin(cell_type)]
sc.pl.dotplot(nk, NK, 'predicted.celltype.l2', dendrogram=False, save="NK")


