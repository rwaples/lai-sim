import pandas as pd
from common.utils import add_reports
from functools import reduce

reports = snakemake.params.reports
out = snakemake.output[0]

reps = [pd.read_csv(x, sep = '\t') for x in reports]
combined_report = reduce(add_reports, reps)
combined_report.to_csv(out, sep ='\t', index=None, float_format='%.3f')
