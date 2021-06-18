import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

site_matrix_path = str(snakemake.input.la_mat)
plot_path = str(snakemake.output.plot)

site_matrix = np.loadz_compressed(file=site_matrix_path)

fig, ax = plt.subplots(figsize = (16,16))
ax = sns.heatmap(site_matrix[::100].T)
plt.title('True local ancestry')
fig.savefig(plot_path)   # save the figure to file
