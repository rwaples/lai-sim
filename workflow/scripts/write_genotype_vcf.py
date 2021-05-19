import tskit
import collections


site_ts = str(snakemake.input.site_ts)
base_path = str(snakemake.params.base_path)
chrom_id = str(WHAT)


# remove the existing 'individuals'  a ts for exporting with inds removed
ts = tszip.decompress(site_ts)
tables = ts.tables
tables.individuals.clear()
tables.nodes.individual = np.repeat(np.array([-1], dtype='int32'), len(tables.nodes))
export_ts = tables.tree_sequence()


## remake a list of individual labels
## labels will work up to 9999 inds per population
pop_of_sample = dict(zip(range(len(export_ts.tables.nodes)), export_ts.tables.nodes.population))
nind_in_pop = collections.defaultdict(int)
pops = [pop_of_sample[i] for i in export_ts.samples()]
for p in pops:
    nind_in_pop[p] +=1
for p in nind_in_pop:
    nind_in_pop[p] = int(nind_in_pop[p]/2)
ind_labels = []
for p in nind_in_pop:
    for i in range(1, nind_in_pop[p]+1):
        ind_labels.append(f'pop_{p}-ind_{i:04}')


with open(base_path + '.genotypes.vcf'), "w") as vcf_file:
	export_ts.write_vcf(
		vcf_file,
		ploidy=2,
		contig_id=chrom_id,
		individual_names=ind_labels,
		position_transform=np.round
	)


## write file mapping inds to populations
with open(base_path + 'sample_map.txt'), 'w') as OUTFILE:
    for ind_string in ind_labels[:nref_samp]:
        pop = ind_string.split('-')[0]
        OUTFILE.write(f'{ind_string}\t{pop}\n')

## write files with inds from the reference and admixed populations


## write a genetic map file
vcf_df = pd.read_csv(base_path + '.genotypes.vcf'), comment='#', sep='\t', header=None)
positions = vcf_df[1].tolist()
rate = contig.recombination_map.get_ll_recombination_map().get_rates()[0] * 100 # convert to cM
cM_pos = [pos*rate for pos in positions]
with open(os.path.join(WORK_DIR, 'genetic_map.txt'), 'w') as OUTFILE:
    for i in range(len(positions)):
        OUTFILE.write(f"{chrom_id}\t{positions[i]}\t{cM_pos[i]:0.6f}\n")


## compress vcf files
!{BCFTOOLS} view -O z \
--samples-file {os.path.join(WORK_DIR, 'admixed_inds.txt')} \
--output {os.path.join(WORK_DIR, "admixed_inds.vcf.gz")} \
{os.path.join(WORK_DIR, "genotypes.bcf")}
## index vcf files
!{BCFTOOLS} index  {base_path + '.genotypes.vcf.gz'}

## delete the raw vcf files
