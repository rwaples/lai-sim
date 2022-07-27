"""Combine via weighted average all QQ reports."""
import pandas as pd
from common.utils import add_reports
from functools import reduce

reports = snakemake.params.reports
out = snakemake.output[0]

reps = []
for x in reports:
	try:
		rep = pd.read_csv(x, sep='\t')
		reps.append(rep)
	except pd.errors.EmptyDataError:
		pass

combined_report = reduce(add_reports, reps)
combined_report.to_csv(out, sep='\t', index=None, float_format='%.4f')
