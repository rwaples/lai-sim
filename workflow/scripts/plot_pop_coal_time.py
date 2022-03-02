import tskit
import tszip
import matplotlib.pyplot as plt
import numpy as np


site_ts = str(snakemake.input.site_ts)
plot_path = str(snakemake.output.plot)
ts = tszip.decompress(site_ts)


for x in range(len(ts.populations())):
	y = ts.tables.nodes.time[np.where(ts.tables.nodes.population == x)[0]]
	plt.plot(np.log10(np.sort(y) + 1), label=x)
plt.legend(title='population')
plt.ylabel('log10(node age+1)')
plt.xlabel('nodes within each population')
plt.savefig(plot_path)
