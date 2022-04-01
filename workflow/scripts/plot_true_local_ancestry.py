import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

site_matrix_path = str(snakemake.input.la_mat)
plot_path = str(snakemake.output.plot)

site_interval = 100  # only plot every 100 sites
nhap = 50  # number of haplotypes to plot

site_matrix = np.load(site_matrix_path)['arr']
fig, ax = plt.subplots(figsize=(32, 16))
ax = sns.heatmap(site_matrix[::site_interval].T[:nhap, :])
plt.title('True local ancestry')
fig.savefig(plot_path)
