import pandas as pd

files = snakemake.params.files
names = snakemake.params.names
nind_ref = snakemake.params.nind_ref
report = snakemake.output[0]

with open (report, 'w') as OUTFILE:
	for i,f in enumerate(files):
		df = pd.read_csv(f, sep ='\t')
		OUTFILE.write(names[i] + '\n')
		OUTFILE.write(nind_ref[i] + '\n')
		df.to_csv(OUTFILE, sep ='\t', index = None)
		OUTFILE.write('\n')
