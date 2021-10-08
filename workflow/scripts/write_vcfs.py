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
target_pop = int(snakemake.params.target_pop)


# remove the existing 'individuals', make a ts for exporting with inds removed
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


# write bcf files with genotypes for all inds
read_fd, write_fd = os.pipe()
write_pipe = os.fdopen(write_fd, "w")
with open(os.path.join(base_path, 'genotypes.bcf'), "w") as bcf_file:
	proc = subprocess.Popen(
	[f"{snakemake.config['PATHS']['BCFTOOLS']}", "view", "-O", "b"], stdin=read_fd, stdout=bcf_file
	)
	export_ts.write_vcf(
		write_pipe,
		ploidy=2,
		contig_id=chrom_id,
		individual_names=ind_labels,
		position_transform=np.round
	)
	write_pipe.close()
	os.close(read_fd)
	proc.wait()
	if proc.returncode != 0:
		raise RuntimeError("bcftools failed with status:", proc.returncode)


# write file mapping inds to populations
# only for the reference individuals
npops = len(export_ts.populations())
nref_total = (npops-1) * nind_ref
with open(os.path.join(base_path, 'sample_map.txt'), 'w') as OUTFILE:
	for ind_string in ind_labels[:nref_total]:
		pop = ind_string.split('-')[0]
		OUTFILE.write(f'{ind_string}\t{pop}\n')

## write files with inds from the reference and admixed populations
with open(os.path.join(base_path, 'reference_inds.txt'), 'w') as OUTFILE:
	for ind_string in ind_labels[:nref_total]:
		OUTFILE.write(f'{ind_string}\n')

with open(os.path.join(base_path, 'target_inds.txt'), 'w') as OUTFILE:
	for ind_string in ind_labels[nref_total:]:
		OUTFILE.write(f'{ind_string}\n')


# write two formats of genetic maps files
os.system(f"{snakemake.config['PATHS']['BCFTOOLS']} query -f '%POS\\n' {os.path.join(base_path, 'genotypes.bcf')} > {os.path.join(base_path, 'site.positions')}")
positions = pd.read_csv(os.path.join(base_path, 'site.positions'), header = None)[0].values
species = stdpopsim.get_species("HomSap")
contig = species.get_contig(chrom_id, length_multiplier=chr_len)
rate = contig.recombination_map.get_rates()[0] * 100 # convert to cM
cM_pos = [pos*rate for pos in positions]
with open(os.path.join(base_path, 'genetic_map.txt') , 'w') as OUTFILE:
	for i in range(len(positions)):
		OUTFILE.write(f"{chrom_id}\t{positions[i]}\t{cM_pos[i]:0.6f}\n")
with open(os.path.join(base_path, 'plink_map.txt'), 'w') as OUTFILE:
	for i in range(len(positions)):
		OUTFILE.write(f"{chrom_id}\t.\t{cM_pos[i]:0.6f}\t{positions[i]}\n")


# write vcfs with true local ancestry
# now this contains only the target indiviudals,  no reference populations
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
header = vcfheader(export_ts, target_pop = target_pop)

read_fd, write_fd = os.pipe()
write_pipe = os.fdopen(write_fd, "w")
with open(os.path.join(base_path, 'la_true.bcf'), "w") as bcf_file:
    proc = subprocess.Popen(
        [f"{snakemake.config['PATHS']['BCFTOOLS']}", "view", "-O", "b"], stdin=read_fd, stdout=bcf_file
    )
    #print(header, file=write_pipe)
    for line in header:
        write_pipe.write(line)

    for i in range(la_mat.shape[0]):
        line = '\t'.join(metadata.iloc[i].astype(str).tolist() + list(np.char.add(np.char.add(la_mat[i, ::2].astype(str),  '|' ), la_mat[i, 1::2].astype(str)))) + '\n'
        write_pipe.write(line)

    write_pipe.close()
    os.close(read_fd)
    proc.wait()
    if proc.returncode != 0:
        raise RuntimeError("bcftools failed with status:", proc.returncode)


## indexing and file conversion
subprocess.run([
	snakemake.config['PATHS']['BCFTOOLS'],
	'index',
	f'{os.path.join(base_path, "genotypes.bcf")}'
	])

subprocess.run([
	snakemake.config['PATHS']['BCFTOOLS'],
	'index',
	f'{os.path.join(base_path, "la_true.bcf")}'
	])

subprocess.run([
	snakemake.config['PATHS']['BCFTOOLS'],
	'view',
	'-O', 'z',
	'--output', f'{os.path.join(base_path, "genotypes.vcf.gz")}',
	f'{os.path.join(base_path, "genotypes.bcf")}'
	])

subprocess.run([
	snakemake.config['PATHS']['BCFTOOLS'],
	'view',
	'-O', 'z',
	'--output', f'{os.path.join(base_path, "la_true.vcf.gz")}',
	f'{os.path.join(base_path, "la_true.bcf")}'
	])

subprocess.run([
	snakemake.config['PATHS']['BCFTOOLS'],
	'index',
	f'{os.path.join(base_path, "genotypes.vcf.gz")}'
	])

subprocess.run([
	snakemake.config['PATHS']['BCFTOOLS'],
	'index',
	f'{os.path.join(base_path, "la_true.vcf.gz")}'
	])

## delete the raw vcf files
#os.remove(f'{os.path.join(base_path, "genotypes.vcf")}')
#os.remove(f'{os.path.join(base_path, "la_true.vcf")}')
