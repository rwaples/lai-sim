import pandas as pd
import numpy as np


rfmix2_fb = str(snakemake.input.inferred_la)
threshold = float(snakemake.params.threshold)
diploid =_path = str(snakemake.output.diploid)

# load rfmix2 forward-backward output
lai = pd.read_csv(rfmix2_fb, comment='#', sep='\t')
lai = lai.set_index(['chromosome',  'genetic_position', 'genetic_marker_index'])
lai = lai.melt(id_vars = ['physical_position'])
lai[['name', 'hap', 'pop']] = lai['variable'].str.split(':::', expand=True)
# apply threshold
lai['calls'] = lai['value'] > threshold

# convert to matrix / dataframe
calls = lai.pivot_table(index = ['physical_position', 'name', 'hap', 'pop'], values = 'calls').reset_index()
lai_calls = calls.query('calls == True').pivot(index = ['physical_position'], columns = ['name', 'hap', 'calls'])
lai_calls = lai_calls.replace('pop_0', 0).replace('pop_1', 1).replace('pop_2', 2).replace('pop_3', 3).replace(np.NaN, -1).astype(int)
del calls
# convert to diploid calls
haploid_calls_1 = lai_calls.values[:, ::2]
haploid_calls_2 = lai_calls.values[:, 1::2]
highpop = np.maximum(lai_calls.values[:, ::2], lai_calls.values[:, 1::2])
lowpop = np.minimum(lai_calls.values[:, ::2], lai_calls.values[:, 1::2])
lai_calls_diploid = lowpop + highpop*10
del haploid_calls_1, haploid_calls_2, highpop, lowpop

# export diploid la calls as data frame
df = pd.DataFrame(lai_calls_diploid)
df.index = lai_calls.index
df.columns = lai_calls.columns.get_level_values(1).values[::2]
df.to_hdf(diploid,
	key = 'la',
	mode = 'w',
	complib = 'blosc:lz4hc',
	complevel = 8)
