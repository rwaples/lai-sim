import pandas as pd

files = snakemake.input
files = [str(x) for x in files]
report = str(snakemake.output)

with open(report, 'w') as OUTFILE:
	for i, f in enumerate(files):
		with open(f) as INFILE:
			for line in INFILE:
				OUTFILE.write(line)
