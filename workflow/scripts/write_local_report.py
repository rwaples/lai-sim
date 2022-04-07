import pandas as pd

files = snakemake.input
files = [str(x) for x in files]
report = str(snakemake.output)

with open(report, 'w') as OUTFILE:
	for i, f in enumerate(files):
		try:
			df = pd.read_csv(f, sep='\t', header=None)
		except pd.errors.EmptyDataError:
			df = pd.DataFrame()
		OUTFILE.write(files[i] + '\t')
		df.to_csv(OUTFILE, sep='\t', index=None, header=None)
