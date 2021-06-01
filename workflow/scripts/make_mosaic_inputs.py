import numpy as np
import pandas as pd
import tszip
import os
from common.utils import make_ind_labels


site_ts = str(snakemake.input.site_ts)
plink_map = str(snakemake.input.plink_map)

folder = str(snakemake.params.folder)
chrom_id = str(snakemake.params.chrom_id)
nind_ref = int(snakemake.params.nind_ref)




# write sample.names
ts = tszip.decompress(site_ts)
npops = len(ts.populations())
ind_labels = make_ind_labels(ts)
nref_total = (npops-1) * nind_ref

with open(os.path.join(folder, 'sample.names'), 'w') as OUTFILE:
	for ind_string in ind_labels[:nref_total]:
		pop = ind_string.split('-')[0]
		OUTFILE.write(f'{pop}\t{ind_string}\n')
	for ind_string in ind_labels[:nref_total]:
		pop = 'admixed'
		OUTFILE.write(f'{pop}\t{ind_string}\n')


# write snpfile
gmap = pd.read_csv(plink_map, sep ='\t', header=None)
gmap.columns = ['chr', 'rsID', 'cM', 'bp']
snpfile = gmap[['rsID' , 'chr', 'cM', 'bp']].copy()
snpfile['A1'] = 'A'
snpfile['A2'] = 'T'
snpfile.to_csv(os.path.join(folder, f'snpfile.{chrom_id}'), sep ='\t', index = None, header=None)


# write rates
with open(os.path.join(folder, f'rates.{chrom_id}'), 'w') as OUTFILE:
	OUTFILE.write(f':sites:{len(snpfile)}\n')
	bp = ' '.join([str(x) for x in snpfile['bp'].values])
	OUTFILE.write(f'{bp}\n')
	cM = ' '.join([str(x) for x in snpfile['cM'].values])
	OUTFILE.write(f'{cM}\n')
