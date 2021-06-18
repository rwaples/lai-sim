import pandas as pd
import numpy as np


site_matrix_path = str(snakemake.input.la_mat)
rfmix2_fb = str(snakemake.input.inferred_la)
threshold = float(snakemake.params.threshold)
stats_path = str(snakemake.output.stats)
plot_path = str(snakemake.output.plot)


# load true local ancestry
loaded = np.load(site_matrix_path)
site_matrix = loaded['arr']

# local rfmix2 forward-backward output
lai = pd.read_csv(rfmix2_fb, comment='#', sep='\t')
lai = lai.set_index(['chromosome',  'genetic_position', 'genetic_marker_index'])
lai = lai.melt(id_vars = ['physical_position'])
lai[['name', 'hap', 'pop']] = lai['variable'].str.split(':::', expand=True)
lai['calls'] = lai['value'] > threshold

calls = lai.pivot_table(index = ['physical_position', 'name', 'hap', 'pop'], values = 'calls').reset_index()
lai_calls = calls.query('calls == True').pivot(index = ['physical_position'], columns = ['name', 'hap', 'calls'])
lai_calls = lai_calls.replace('pop_0', 0).replace('pop_1', 1).replace('pop_2', 2).replace('pop_3', 3).replace(np.NaN, -1).astype(int)
del calls
haploid_calls_1 = lai_calls.values[:, ::2]
haploid_calls_2 = lai_calls.values[:, 1::2]
highpop = np.maximum(lai_calls.values[:, ::2], lai_calls.values[:, 1::2])
lowpop = np.minimum(lai_calls.values[:, ::2], lai_calls.values[:, 1::2])
lai_calls_diploid = lowpop + highpop*10
del haploid_calls_1, haploid_calls_2, highpop, lowpop
