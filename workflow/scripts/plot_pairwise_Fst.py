import tskit
import tszip
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


site_ts = str(snakemake.input.site_ts)
plot_path = str(snakemake.output.plot)
ts = tszip.decompress(site_ts)

npops = len(ts.populations())


def Hudson_Fst(pop1, pop2, ts):
	pop1_samples = ts.samples(population=pop1)
	pop2_samples = ts.samples(population=pop2)
	# once the two new sets of samples are specified, proceed as above
	d1 = ts.diversity(pop1_samples, windows = 'sites', span_normalise=False)
	d2 = ts.diversity(pop2_samples, windows = 'sites', span_normalise=False)
	d12 = ts.divergence([pop1_samples, pop2_samples], windows = 'sites', span_normalise=False)
	mean_within = (d1 + d2)/2
	Fst = 1 - mean_within.sum()/d12.sum()
	return(Fst)

mat = np.zeros((npops, npops))
mat.fill(np.NaN)
for pop1 in range(npops):
	for pop2 in range(pop1+1, npops):
		Fst = Hudson_Fst(pop1, pop2, ts)
		mat[pop1, pop2] = Fst
		mat[pop2, pop1] = Fst

plot = sns.heatmap(mat, annot=True, annot_kws={'fontsize':12})
plt.title('Fst')
plt.savefig(plot_path)
