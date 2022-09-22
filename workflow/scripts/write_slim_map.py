"""Write a genetic map file for use by SLiM."""
import msprime
import numpy as np
import pandas as pd

rec_map_path = str(snakemake.params.rec_map_path)
slim_map_path = str(snakemake.output.slim_map_path)
#max_bp = int(snakemake.params.max_bp)

rmap = msprime.RateMap.read_hapmap(rec_map_path)

map_df = pd.DataFrame([rmap.left, rmap.right, rmap.span, rmap.rate]).T
map_df.columns = ['left', 'right', 'span', 'rate']

# fill any missing values with the mean rate
mean_rate = rmap.mean_rate
map_df = map_df.fillna(rmap.mean_rate)

# cut off the map if it extends past @max_bp
# also cut of the span from 0 to 1
#if map_df['right'].max() < max_bp:
#	shorted = map_df.query("right>1").query("left<@max_bp").copy()
#	ends = shorted['right'].values.astype(int)
#	ends[-1] = max_bp
#	rates = shorted['rate'].values
#else:

map_df = map_df.query("right>1").copy()
ends = map_df['right'].values.astype(int)
rates = map_df['rate'].values

assert len(rates) == len(ends)

# adjust to zero-based indexing
ends = ends - 1

with open(slim_map_path, 'w') as OUTFILE:
	for i in range(len(ends)):
		OUTFILE.write(f'{ends[i]}\t{rates[i]}\n')
