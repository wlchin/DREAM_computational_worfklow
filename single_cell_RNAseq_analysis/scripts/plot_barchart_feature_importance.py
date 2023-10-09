import pandas as pd
import matplotlib.pyplot as plt

ranks_df = pd.read_csv("data/top_rank_loadings_df_blood.csv")
R_genes = pd.read_csv("results/blood_markers_for_dotplot_R.txt", header = None)
NR_genes = pd.read_csv("results/blood_markers_for_dotplot_NR.txt", header = None)

res_only = ranks_df[ranks_df["Unnamed: 0"].isin(R_genes[0])]
non_res_only = ranks_df[ranks_df["Unnamed: 0"].isin(NR_genes[0])]

combined_df = pd.concat([res_only, non_res_only], axis = 0)
sorted_df = combined_df.sort_values("value.var")

categories = sorted_df["Unnamed: 0"].to_list()
values = sorted_df["value.var"] * -1

data = list(zip(categories, values))
sorted_data = sorted(data, key=lambda x: x[1])
sorted_categories, sorted_values = zip(*sorted_data)
# Create a figure and axis
fig, ax = plt.subplots(figsize = (2,12.5))

# Create the barplot
bars = ax.barh(sorted_categories, sorted_values, color=['blue' if v >= 0 else 'red' for v in values])

# Add labels and title
ax.set_xlabel('Loading scores')
ax.set_ylabel('Genes')
ax.set_title('Feature importance')
# Add value labels above the bars

# Show the plot

plt.savefig('results/barplot_feature_importance.png', bbox_inches='tight', dpi = 600)