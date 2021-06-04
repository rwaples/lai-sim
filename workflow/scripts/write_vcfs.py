import tskit
import tszip
import stdpopsim
import numpy as np
import pandas as pd
import collections
import subprocess
import os
from common.utils import make_ind_labels, vcfheader


site_ts = str(snakemake.input.site_ts)
site_matrix = str(snakemake.input.site_matrix)

base_path = str(snakemake.params.base_path)
chrom_id = str(snakemake.params.chrom_id)
nind_ref = int(snakemake.params.nind_ref)
chr_len = float(snakemake.params.chr_len)

# remove the existing 'individuals'  a ts for exporting with inds removed
ts = tszip.decompress(site_ts)
tables = ts.tables
tables.individuals.clear()
tables.nodes.individual = np.repeat(np.array([-1], dtype='int32'), len(tables.nodes))
# set ancestral alleles to 'A', derived alleles to 'T'
a, off = tskit.pack_strings(['A' for _ in tables.sites])
tables.sites.set_columns(
	position=tables.sites.position,
	ancestral_state=a,
	ancestral_state_offset=off
	)
t, off = tskit.pack_strings(['T' for _ in tables.mutations])
tables.mutations.set_columns(
	site=tables.mutations.site,
	node=tables.mutations.node,
	derived_state=t,
	derived_state_offset=off
	)
export_ts = tables.tree_sequence()

ind_labels = make_ind_labels(export_ts)

with open(base_path + '.genotypes.vcf', "w") as vcf_file:
	export_ts.write_vcf(
		vcf_file,
		ploidy=2,
		contig_id=chrom_id,
		individual_names=ind_labels,
		position_transform=np.round
	)

npops = len(export_ts.populations())
nref_total = (npops-1) * nind_ref


## write file mapping inds to populations
# only for the reference individuals
with open(base_path + '.sample_map.txt', 'w') as OUTFILE:
    for ind_string in ind_labels[:nref_total]:
        pop = ind_string.split('-')[0]
        OUTFILE.write(f'{ind_string}\t{pop}\n')

## write files with inds from the reference and admixed populations
with open(base_path + '.reference_inds.txt', 'w') as OUTFILE:
    for ind_string in ind_labels[:nref_total]:
        #pop = ind_string.split('-')[0]
        OUTFILE.write(f'{ind_string}\n')

with open(base_path + '.target_inds.txt', 'w') as OUTFILE:
    for ind_string in ind_labels[nref_total:]:
        #pop = ind_string.split('-')[0]
        OUTFILE.write(f'{ind_string}\n')


## write a genetic map file
species = stdpopsim.get_species("HomSap")
contig = species.get_contig(chrom_id, length_multiplier=chr_len)
rate = contig.recombination_map.get_rates()[0] * 100 # convert to cM
vcf_df = pd.read_csv(base_path + '.genotypes.vcf', comment='#', sep='\t', header=None)
positions = vcf_df[1].tolist()
del vcf_df
cM_pos = [pos*rate for pos in positions]
with open(base_path + '.genetic_map.txt', 'w') as OUTFILE:
    for i in range(len(positions)):
        OUTFILE.write(f"{chrom_id}\t{positions[i]}\t{cM_pos[i]:0.6f}\n")

with open(base_path + '.plink_map.txt', 'w') as OUTFILE:
    for i in range(len(positions)):
        OUTFILE.write(f"{chrom_id}\t.\t{cM_pos[i]:0.6f}\t{positions[i]}\n")


## write vcfs with true local ancestry
la_mat = np.load(file = site_matrix)['arr']
df = pd.DataFrame(la_mat)
df.insert(loc = 0, column= 'CHROM', value = chrom_id)
df.insert(loc = 1, column= 'POS', value = positions)
df.insert(loc = 2, column= 'ID', value = [f'site_{x.id}' for x in export_ts.sites()])
df.insert(loc = 3, column= 'REF', value = 'A')
df.insert(loc = 4, column= 'ALT', value = 'T')
df.insert(loc = 5, column= 'QUAL', value = '.')
df.insert(loc = 6, column= 'FILTER', value = 'PASS')
df.insert(loc = 7, column= 'INFO', value = '.')
df.insert(loc = 8, column= 'FORMAT', value = 'LA')
metadata = df.iloc[:,0:9].copy()
for i in range(len(ind_labels)):
	diploid_anc = df.loc[:,i*2].astype(str) + '|' + df.loc[:,i*2+1].astype(str)
	metadata[ind_labels[i]] = diploid_anc

header = vcfheader(
	contig_id=chrom_id,
	contig_len=int(np.ceil(export_ts.get_sequence_length())),
	ind_labels=ind_labels
	)

output_VCF = base_path + ".la_true.vcf"
with open(output_VCF, 'w') as vcf:
    vcf.write(header)

metadata.to_csv(output_VCF, sep="\t", mode='a', index=False, header=False)


## use BCFTOOLS to compress vcf files with gz
subprocess.run([
	snakemake.config['PATHS']['BCFTOOLS'],
	'view',
	'-O', 'z',
	'--output', f'{base_path}.genotypes.vcf.gz',
	f'{base_path}.genotypes.vcf'
	])

subprocess.run([
	snakemake.config['PATHS']['BCFTOOLS'],
	'view',
	'-O', 'z',
	'--output', f'{base_path}.la_true.vcf.gz',
	f'{base_path}.la_true.vcf'
	])


subprocess.run([
	snakemake.config['PATHS']['BCFTOOLS'],
	'index',
	f'{base_path}.genotypes.vcf.gz'
	])

subprocess.run([
	snakemake.config['PATHS']['BCFTOOLS'],
	'index',
	f'{base_path}.la_true.vcf.gz'
	])

## delete the raw vcf files
os.remove(f'{base_path}.genotypes.vcf')
os.remove(f'{base_path}.la_true.vcf')
