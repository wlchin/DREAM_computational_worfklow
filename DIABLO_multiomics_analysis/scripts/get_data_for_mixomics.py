
import pandas as pd

number_of_genes = 10000

nanostring = pd.read_csv(snakemake.input[0], index_col = 0)
cols = []

for i in nanostring.columns:
    if len(i) == 4:
        num = i[2:4]
        new_label = "PT" + num
        cols.append(new_label)
        
    if len(i) == 3:
        num = i[2:]
        new_label = "PT0" + num
        cols.append(new_label)

nanostring.columns = cols
nanostring_all = nanostring.T

pheno = pd.read_csv(snakemake.input[1])
pheno_small = pheno[["patient", "pd_6m"]].drop_duplicates()
pheno_small.index = pheno_small.patient

totalmat = pd.read_csv(snakemake.input[2], index_col = 0)

baseline_mat = totalmat[totalmat.columns[totalmat.columns.str.contains("Baseline")]].T
baseline_mat.index = [i[:4] for i in baseline_mat.index]

commons = set(baseline_mat.index).intersection(set(nanostring_all.index)) # 34 samples common to both
#pheno_small["pd_6m"][pheno_small.index.isin(commons)].value_counts() # not bad representation - check balance

phenodata = pheno_small[pheno_small.index.isin(commons)]
phenodata.to_csv("data/phenodata_PLS.csv")

tumour = nanostring_all[nanostring_all.index.isin(commons)].reindex(phenodata.index)
tumour.to_csv("data/tumour_PLS.csv")

baseline_mat = totalmat[totalmat.columns[totalmat.columns.str.contains("Baseline")]].T
baseline_mat.index = [i[:4] for i in baseline_mat.index]
peripheral_blood = baseline_mat[baseline_mat.index.isin(commons)].reindex(phenodata.index)
topgenes = peripheral_blood[peripheral_blood.columns[peripheral_blood.mean(0) > 1]].var(0).sort_values().tail(number_of_genes).index
peripheral_blood[topgenes].to_csv("data/peripheralbloodfiltered_PLS.csv")

C2D1_mat = totalmat[totalmat.columns[totalmat.columns.str.contains("C2D1")]].T
C2D1_mat.index = [i[:4] for i in baseline_mat.index]
peripheral_blood = C2D1_mat[C2D1_mat.index.isin(commons)].reindex(phenodata.index)
topgenes = peripheral_blood[peripheral_blood.columns[peripheral_blood.mean(0) > 1]].var(0).sort_values().tail(number_of_genes).index
peripheral_blood[topgenes].to_csv("data/peripheralbloodC2D1filtered_PLS.csv")

baseline_mat = totalmat[totalmat.columns[totalmat.columns.str.contains("C3D1")]].T
baseline_mat.index = [i[:4] for i in baseline_mat.index]
peripheral_blood = baseline_mat[baseline_mat.index.isin(commons)].reindex(phenodata.index)
topgenes = peripheral_blood[peripheral_blood.columns[peripheral_blood.mean(0) > 1]].var(0).sort_values().tail(number_of_genes).index
peripheral_blood[topgenes].to_csv("data/peripheralbloodC3D1filtered_PLS.csv")