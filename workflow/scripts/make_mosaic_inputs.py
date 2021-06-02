import numpy as np
import pandas as pd
import tszip
import os
import subprocess
from common.utils import make_ind_labels


target_phased_vcf = str(snakemake.input.target_phased_vcf)
reference_vcf = str(snakemake.input.reference_vcf)

site_ts = str(snakemake.input.site_ts)
plink_map = str(snakemake.input.plink_map)

folder = str(snakemake.params.folder)
chrom_id = str(snakemake.params.chrom_id).strip('chr')
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
	for ind_string in ind_labels[nref_total:]:
		pop = 'admixed'
		OUTFILE.write(f'{pop}\t{ind_string}\n')

# needed to subset vcf files
for p in range(npops-1):
	pop = 'pop_' + (str(p))
	with open(os.path.join(folder, f'{pop}.samplelist'), 'w') as OUTFILE:
		for ind_string in ind_labels[nind_ref*p:nind_ref*(p+1)]:
			OUTFILE.write(f'{ind_string}\n')
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
		'--output', f'{output}',
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
	'--output', f'{output}',
	f'{target_phased_vcf}'
	])

subprocess.run([
	snakemake.config['PATHS']['BCFTOOLS'],
	'index', f'{output}',
	])

# write genofiles
for p in range(npops-1):
	pop = 'pop_' + (str(p))
	prefix = os.path.join(folder, f'{pop}.genofile.{chrom_id}')
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
prefix = os.path.join(folder, f'{pop}.genofile.{chrom_id}')
subprocess.run([
	snakemake.config['PATHS']['BCFTOOLS'],
	'convert',
	'--haplegendsample', f'{prefix}',
	f'{output}',
	])
print(f"zcat {prefix+'.hap.gz'} | tr -d ' ' > {prefix}")
subprocess.run(
	f"zcat {prefix+'.hap.gz'} | tr -d ' ' > {prefix}", shell=True)





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
