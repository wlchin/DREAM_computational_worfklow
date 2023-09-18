import scanpy as sc

blood = sc.read(snakemake.input[0])
blood.var_names = blood.var["gene"]

sc.pp.filter_cells(blood, min_genes=200)
blood.var['mt'] = blood.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(blood, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
blood = blood[blood.obs.pct_counts_mt < 10, :]

sc.pp.normalize_total(blood, target_sum=1e4)
sc.pp.log1p(blood)
blood.raw = blood

sc.pp.highly_variable_genes(blood, n_top_genes = 2500, min_mean=0.0125, max_mean=3, min_disp=0.5)
blood2 = blood[: , blood.var.highly_variable]
sc.pp.scale(blood2, max_value=10)
sc.tl.pca(blood2, svd_solver='arpack')
sc.pp.neighbors(blood2)
sc.tl.umap(blood2)
sc.tl.leiden(blood2)

blood2.write(snakemake.output[0])
blood2.obs.to_csv(snakemake.output[1])