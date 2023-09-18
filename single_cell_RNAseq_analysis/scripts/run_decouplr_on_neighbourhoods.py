import scanpy as sc
import decoupler as dc
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import pandas as pd
import numpy as np
from sklearn.feature_selection import VarianceThreshold

adata_neighbourhoods = sc.read(snakemake.input[0])
msigdb = pd.read_pickle(snakemake.input[1])

dc.run_ora(mat=adata_neighbourhoods, net=msigdb, source='geneset', target='genesymbol', verbose=True, use_raw=False)

acts = dc.get_acts(adata_neighbourhoods, obsm_key='ora_estimate')

hallmarks = acts.obsm["ora_estimate"]

sel = VarianceThreshold(threshold=8)

acts_v = hallmarks.to_numpy()
max_e = np.nanmax(acts_v[np.isfinite(acts_v)])
hallmarks[~np.isfinite(hallmarks)] = max_e
hallmarks_pathways = hallmarks.transpose()
sel.fit_transform(hallmarks)
selected_hallmarks = sel.get_feature_names_out()

hallmarks[selected_hallmarks].to_csv(snakemake.output[0])