"""Write the genetic map to files in two formats used by downstream analysis programs."""
import msprime
import pandas as pd
import numpy as np

sites_path = str(snakemake.input.sites_file)
rec_map_path = str(snakemake.params.rec_map_path)
genetic_map_path = str(snakemake.output.genetic_map)
plink_map_path = str(snakemake.output.plink_map)

chrom_id = '22'

# write two formats of genetic maps files
positions = pd.read_csv(sites_path, header=None)[0].values
rmap = msprime.RateMap.read_hapmap(rec_map_path)

# need to load in the recombination rate map here.
# for the positions in the analyis, I need to be able to calculate the cM position.

# factor of 100 is to convert Morgans -> centiMorgans
cM_pos = np.array([rmap.get_cumulative_mass(x) * 100 for x in positions])

with open(genetic_map_path, 'w') as OUTFILE:
	for i in range(len(positions)):
		OUTFILE.write(f'chr{chrom_id}\t{positions[i]}\t{cM_pos[i]:0.6f}\n')

with open(plink_map_path, 'w') as OUTFILE:
	for i in range(len(positions)):
		OUTFILE.write(f'chr{chrom_id}\t.\t{cM_pos[i]:0.6f}\t{positions[i]}\n')
