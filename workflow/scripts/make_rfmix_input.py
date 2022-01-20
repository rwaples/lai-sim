import tskit
import tszip
import numpy as np
import pandas as pd

site_ts = str(snakemake.input.site_ts)
genetic_map = str(snakemake.input.genetic_map)

alleles_file = str(snakemake.output.alleles_file)
classes_file = str(snakemake.output.classes_file)
positions_file = str(snakemake.output.positions_file)

nind_admixed = int(snakemake.params.nind_admixed)
nind_ref = str(snakemake.params.nind_ref)

nind_ref = np.array([int(x) for x in nind_ref.split(',')])

# write alleles file - ne row per site, one column per haplotype
ts = tszip.decompress(site_ts)
np.savetxt(
	fname=alleles_file,
	X=ts.genotype_matrix(),
	delimiter='',
	fmt = '%i1'
)

# write classes file
npops = len(ts.populations())
classes = list(np.repeat(np.arange(1, npops), 2*nind_ref)) + [0]*2*nind_admixed
classes = [str(x) for x in classes]
with open(classes_file, 'w') as OUTFILE:
	OUTFILE.write(' '.join(classes))
	OUTFILE.write('\n')

# write positions_file - 1 row per site - pos in cM
positions = pd.read_csv(genetic_map, sep = '\t', header = None)[2].values
with open(positions_file, 'w') as OUTFILE:
	for p in positions:
		OUTFILE.write(f"{p:0.6f}\n")
