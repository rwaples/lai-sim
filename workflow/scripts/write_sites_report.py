import pandas as pd

files = snakemake.input
files = [str(x) for x in files]
report = str(snakemake.output)

with open(report, 'w') as OUTFILE:
	for f in files:
		with open(f, 'r') as fp:
			for count, _ in enumerate(fp):
				pass
		nsites = count + 1
		OUTFILE.write(f'{f}\t{nsites}\n')
