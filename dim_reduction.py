from chemplot import Plotter
from data_fetcher import df

# Chemplot
# =========
cp = Plotter.from_smiles(df.canonical_smiles, target=df.max_epo_conc_mean_6hr, sim_type='structural')
df_pca = cp.pca()
custom_plot_pca = cp.interactive_plot(kind='scatter')

cp = Plotter.from_smiles(df.canonical_smiles, target=df.max_epo_conc_mean_6hr, sim_type='structural')
df_umap = cp.umap()
custom_plot_umap = cp.interactive_plot(kind='scatter')

#df_2_components = cp.tsne(random_state=random_state)