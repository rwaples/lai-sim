import numpy as np
import pandas as pd
import tszip
import os
import subprocess
import csv
from common.utils import make_ind_labels


target_phased_vcf = str(snakemake.input.target_phased_vcf)
reference_vcf = str(snakemake.input.reference_vcf)

site_ts = str(snakemake.input.site_ts)
plink_map = str(snakemake.input.plink_map)

folder = str(snakemake.params.folder)
chrom_id = str(snakemake.params.chrom_id).strip('chr')
nind_ref = str(snakemake.params.nind_ref)

nind_ref = np.array([int(x) for x in nind_ref.split(',')])

# write sample.names
ts = tszip.decompress(site_ts)
npops = len(ts.populations())
ind_labels = make_ind_labels(ts)
nref_total = nind_ref.sum()
with open(os.path.join(folder, 'sample.names'), 'w') as OUTFILE:
	for ind_string in ind_labels[:nref_total]:
		pop = ind_string.split('-')[0]
		OUTFILE.write(f'{pop} {ind_string} 0 0 0 2 -9\n')
	for ind_string in ind_labels[nref_total:]:
		pop = 'admixed'
		OUTFILE.write(f'{pop} {ind_string} 0 0 0 2 -9\n')

# needed to subset vcf files
start_ix = 0
stop_ix = 0
for p in range(npops-1):
	pop = 'pop_' + (str(p))
	stop_ix += nind_ref[p]
	with open(os.path.join(folder, f'{pop}.samplelist'), 'w') as OUTFILE:
		for ind_string in ind_labels[start_ix:stop_ix]:
			OUTFILE.write(f'{ind_string}\n')
	start_ix += nind_ref[p]

with open(os.path.join(folder, 'admixed.samplelist'), 'w') as OUTFILE:
	for ind_string in ind_labels[nref_total:]:
		OUTFILE.write(f'{ind_string}\n')

for p in range(npops-1):
	pop = 'pop_' + (str(p))
	samplelist = os.path.join(folder, f'{pop}.samplelist')
	output = os.path.join(folder, f'{pop}.phased.vcf.gz')

	subprocess.run([
		snakemake.config['PATHS']['BCFTOOLS'],
		'view',
		'--samples-file', f'{samplelist}',
		'--output-type', 'z',
		'--output-file', f'{output}',
		f'{reference_vcf}'
		])
	subprocess.run([
		snakemake.config['PATHS']['BCFTOOLS'],
		'index', f'{output}',
		])


pop = 'admixed'
samplelist = os.path.join(folder, f'{pop}.samplelist')
output = os.path.join(folder, f'{pop}.phased.vcf.gz')
subprocess.run([
	snakemake.config['PATHS']['BCFTOOLS'],
	'view',
	'--samples-file', f'{samplelist}',
	'--output-type', 'z',
	'--output-file', f'{output}',
	f'{target_phased_vcf}'
	])

subprocess.run([
	snakemake.config['PATHS']['BCFTOOLS'],
	'index', f'{output}',
	])

# write genofiles
for p in range(npops-1):
	pop = 'pop_' + (str(p))
	prefix = os.path.join(folder, f'{pop}genofile.{chrom_id}')
	output = os.path.join(folder, f'{pop}.phased.vcf.gz')
	subprocess.run([
		snakemake.config['PATHS']['BCFTOOLS'],
		'convert',
		'--haplegendsample', f'{prefix}',
		f'{output}',
		])
	print(f"zcat {prefix+'.hap.gz'} | tr -d ' ' > {prefix}")
	subprocess.run(
		f"zcat {prefix+'.hap.gz'} | tr -d ' ' > {prefix}", shell=True)


pop = 'admixed'
prefix = os.path.join(folder, f'{pop}genofile.{chrom_id}')
output = os.path.join(folder, f'{pop}.phased.vcf.gz')
subprocess.run([
	snakemake.config['PATHS']['BCFTOOLS'],
	'convert',
	'--haplegendsample', f'{prefix}',
	f'{output}',
	])
#print(f"zcat {prefix+'.hap.gz'} | tr -d ' ' > {prefix}")
subprocess.run(
	f"zcat {prefix+'.hap.gz'} | tr -d ' ' > {prefix}", shell=True)


# write snpfile
gmap = pd.read_csv(plink_map, sep ='\t', header=None)
gmap.columns = ['chr', 'rsID', 'cM', 'bp']
snpfile = gmap[['rsID' , 'chr', 'cM', 'bp']].copy()
snpfile['rsID'] = [f'        mosaic_SNP_{i}' for i in range(len(snpfile))]
snpfile['chr'] = [f'{c[3:]}' for c in snpfile['chr']]
snpfile['A1'] = 'A'
snpfile['A2'] = 'T'
snpfile.to_csv(os.path.join(folder, f'snpfile.{chrom_id}'), sep =' ',
	index=None, header=None, quoting=csv.QUOTE_NONE, escapechar = ' ')


# write rates
with open(os.path.join(folder, f'rates.{chrom_id}'), 'w') as OUTFILE:
	OUTFILE.write(f':sites:{len(snpfile)}\n')
	bp = ' '.join([str(x) for x in snpfile['bp'].values])
	OUTFILE.write(f'{bp}\n')
	cM = ' '.join([str(x) for x in snpfile['cM'].values])
	OUTFILE.write(f'{cM}\n')
