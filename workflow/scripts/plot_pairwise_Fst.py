"""Plot a heatmap of pairwise Fst (Hudson)."""
import tskit
import tszip
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from common.utils import get_fst_faser

site_ts = str(snakemake.input.site_ts)
plot_path = str(snakemake.output.plot)
ts = tszip.decompress(site_ts)

npops = len(ts.populations())

mat = np.zeros((npops, npops))
mat.fill(np.NaN)
for pop1 in range(npops):
	for pop2 in range(pop1 + 1, npops):
		Fst = get_fst_faser(ts, pop1, pop2)
		mat[pop1, pop2] = Fst
		mat[pop2, pop1] = Fst

plot = sns.heatmap(mat, annot=True, annot_kws={'fontsize': 12})
plt.title('Fst')
plt.savefig(plot_path)
