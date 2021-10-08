rule la_and_vcf_done:
	input:
		[f'results/{u.model_name}/{u.sim_name}/{u.anal_name}/true_local_ancestry.tracts.hdf' for u in units.itertuples()],
		[f'results/{u.model_name}/{u.sim_name}/{u.anal_name}/true_local_ancestry.site_matrix.npz' for u in units.itertuples()],
		[f'results/{u.model_name}/{u.sim_name}/{u.anal_name}/true_local_ancestry.samples.txt' for u in units.itertuples()],
		[f'results/{u.model_name}/{u.sim_name}/{u.anal_name}/la_true.vcf.gz' for u in units.itertuples()],
		[f'results/{u.model_name}/{u.sim_name}/{u.anal_name}/genotypes.witherror.vcf.gz' for u in units.itertuples()],


rule write_vcfs:
	input:
		site_matrix = 'results/{model_name}/{sim_name}/{anal_name}/true_local_ancestry.site_matrix.npz',
		site_ts = 'results/{model_name}/{sim_name}/{anal_name}/sample.filter.tsz',
	output:
		phased_ancestry_vcf = 'results/{model_name}/{sim_name}/{anal_name}/la_true.vcf.gz',
		phased_genotype_vcf = 'results/{model_name}/{sim_name}/{anal_name}/genotypes.vcf.gz',
		genetic_map = 'results/{model_name}/{sim_name}/{anal_name}/genetic_map.txt',
		plink_map = 'results/{model_name}/{sim_name}/{anal_name}/plink_map.txt',
		sample_map = 'results/{model_name}/{sim_name}/{anal_name}/sample_map.txt',
		reference_inds = 'results/{model_name}/{sim_name}/{anal_name}/reference_inds.txt',
		target_inds = 'results/{model_name}/{sim_name}/{anal_name}/target_inds.txt',
	params:
		base_path = 'results/{model_name}/{sim_name}/{anal_name}',
		chrom_id = lambda w: units.loc[(w.sim_name, w.anal_name)].chr,
		nind_ref = lambda w: units.loc[(w.sim_name, w.anal_name)].nind_ref,
		chr_len = lambda w: units.loc[w.sim_name, w.anal_name].chr_len,
		target_pop = lambda w: units.loc[w.sim_name, w.anal_name].target_pop,
	script:
		'../scripts/write_vcfs.py'


rule split_vcf:
	# split vcfs after phasing
	input:
		vcf = 'results/{model_name}/{sim_name}/{anal_name}/phased.vcf.gz',
		samples = 'results/{model_name}/{sim_name}/{anal_name}/{sample_group}.txt',
	output:
		vcf = 'results/{model_name}/{sim_name}/{anal_name}/phased.{sample_group}.vcf.gz',
	params:
		bcftools = config['PATHS']['BCFTOOLS']
	shell:
		"""
		{params.bcftools} view --samples-file {input.samples} --output {output.vcf} --output-type z {input.vcf}
		{params.bcftools} index {output.vcf}
		"""

rule split_bcf:
	input:
		vcf =  'results/{model_name}/{sim_name}/{anal_name}/phased.vcf.gz',
		samples = 'results/{model_name}/{sim_name}/{anal_name}/{sample_group}.txt'
	output:
		bcf = 'results/{model_name}/{sim_name}/{anal_name}/phased.{sample_group}.bcf',
	params:
		bcftools = config['PATHS']['BCFTOOLS']
	shell:
		"""
		{params.bcftools} view --samples-file {input.samples} --output {output.bcf} --output-type b {input.vcf}
		{params.bcftools} index {output.bcf}
		"""


rule true_local_ancestry:
	input:
		ancestry_ts = 'results/{model_name}/{sim_name}/{anal_name}/sample.tsz',
		site_ts = 'results/{model_name}/{sim_name}/{anal_name}/sample.filter.tsz',
		node_mapping = 'results/{model_name}/{sim_name}/{anal_name}/sample.filter.node_mapping.npy'
	output:
		tracts = 'results/{model_name}/{sim_name}/{anal_name}/true_local_ancestry.tracts.hdf',
		site_matrix = 'results/{model_name}/{sim_name}/{anal_name}/true_local_ancestry.site_matrix.npz',
		samples = 'results/{model_name}/{sim_name}/{anal_name}/true_local_ancestry.samples.txt',
	params:
		admixture_time = lambda w: units.loc[(w.sim_name, w.anal_name)].admixture_time,
		target_pop = lambda w: units.loc[(w.sim_name, w.anal_name)].target_pop,
		nind_admixed = lambda w: units.loc[(w.sim_name, w.anal_name)].nind_admixed,
	script:
		'../scripts/true_local_ancestry.py'


rule sample_sites:
	input:
		'results/{model_name}/{sim_name}/full.tsz'
	output:
		'results/{model_name}/{sim_name}/{anal_name}/sample.tsz'
	params:
		nind_admixed = lambda w: units.loc[(w.sim_name, w.anal_name)].nind_admixed,
		nind_ref = lambda w: units.loc[(w.sim_name, w.anal_name)].nind_ref,
		anal_seed = lambda w: units.loc[(w.sim_name, w.anal_name)].anal_seed,
		admixture_time = lambda w: units.loc[(w.sim_name, w.anal_name)].admixture_time,
	script:
		"../scripts/sample.py"


rule filter_inds:
	input:
		'results/{model_name}/{sim_name}/{anal_name}/sample.tsz'
	output:
		tsz = 'results/{model_name}/{sim_name}/{anal_name}/sample.filter.tsz',
		node_mapping = 'results/{model_name}/{sim_name}/{anal_name}/sample.filter.node_mapping.npy'
	params:
		MAC_filter = lambda w: units.loc[(w.sim_name, w.anal_name)].MAC_filter,
		max_snps = lambda w: units.loc[(w.sim_name, w.anal_name)].max_snps,
		anal_seed = lambda w: units.loc[(w.sim_name, w.anal_name)].anal_seed,
	script:
		"../scripts/filter.py"
