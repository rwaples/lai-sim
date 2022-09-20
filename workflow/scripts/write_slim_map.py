"""Write a genetic map file for use by SLiM."""
import msprime
import numpy as np
import pandas as pd

rec_map_path = str(snakemake.params.rec_map_path)
slim_map_path = str(snakemake.output.slim_map_path)
max_bp = int(snakemake.params.max_bp)

rmap = msprime.RateMap.read_hapmap(rec_map_path)

map_df = pd.DataFrame([rmap.left, rmap.right, rmap.span, rmap.rate]).T
map_df.columns = ['left', 'right', 'span', 'rate']
# fill any missing values with the mean rate
indices = ~np.isnan(map_df['rate'].values)
mean_rate = np.average(
	map_df['rate'].values[indices],
	weights=map_df['span'].values[indices]
)
map_df = map_df.fillna(mean_rate)
shorted = map_df.query("left<@max_bp").copy()
ends = shorted['right'].values.astype(int)
ends[-1] = max_bp
rates = shorted['rate'].values
assert len(rates) == len(ends)
with open(slim_map_path, 'w') as OUTFILE:
	for i in range(len(ends)):
		OUTFILE.write(f'{ends[i]}\t{rates[i]}\n')
