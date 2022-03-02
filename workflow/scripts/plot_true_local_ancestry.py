import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

site_matrix_path = str(snakemake.input.la_mat)
plot_path = str(snakemake.output.plot)

loaded = np.load(site_matrix_path)
site_matrix = loaded['arr']
fig, ax = plt.subplots(figsize=(32, 16))
ax = sns.heatmap(site_matrix[::100].T)  # only plot every 100 sites
plt.title('True local ancestry')
fig.savefig(plot_path)
