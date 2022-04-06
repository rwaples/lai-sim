import pandas as pd

files = snakemake.params.input
report = snakemake.output

with open(report, 'w') as OUTFILE:
	for i, f in enumerate(files):
		# df = pd.read_csv(f, sep='\t')
		# OUTFILE.write(names[i] + '\n')
		# OUTFILE.write(nind_ref[i] + '\n')
		OUTFILE.write(files[i] + '\t')
		df.to_csv(OUTFILE, sep='\t', index=None)
		OUTFILE.write('\n')
