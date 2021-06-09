rule la_and_vcf_done:
	input:
		[f'results/local_ancestry/{u.model_name}/{u.sim_name}/{u.anal_name}.true_local_ancestry.tracts.hdf' for u in units.itertuples()],
		[f'results/local_ancestry/{u.model_name}/{u.sim_name}/{u.anal_name}.true_local_ancestry.site_matrix.npz' for u in units.itertuples()],
		[f'results/local_ancestry/{u.model_name}/{u.sim_name}/{u.anal_name}.true_local_ancestry.samples.txt' for u in units.itertuples()],
		[f'results/local_ancestry/{u.model_name}/{u.sim_name}/{u.anal_name}.la_true.vcf.gz' for u in units.itertuples()],
		[f'results/local_ancestry/{u.model_name}/{u.sim_name}/{u.anal_name}.genotypes.vcf.gz' for u in units.itertuples()],


rule write_vcfs:
	input:
		site_matrix = 'results/local_ancestry/{model_name}/{sim_name}/{anal_name}.true_local_ancestry.site_matrix.npz',
		site_ts = 'results/simulations/{model_name}/{sim_name}/{anal_name}.sample.filter.tsz',
	output:
		phased_ancestry_vcf = 'results/local_ancestry/{model_name}/{sim_name}/{anal_name}.la_true.vcf.gz',
		phased_genotype_vcf = 'results/local_ancestry/{model_name}/{sim_name}/{anal_name}.genotypes.vcf.gz',
		genetic_map = 'results/local_ancestry/{model_name}/{sim_name}/{anal_name}.genetic_map.txt',
		plink_map = 'results/local_ancestry/{model_name}/{sim_name}/{anal_name}.plink_map.txt',
		sample_map = 'results/local_ancestry/{model_name}/{sim_name}/{anal_name}.sample_map.txt',
		reference_inds = 'results/local_ancestry/{model_name}/{sim_name}/{anal_name}.reference_inds.txt',
		target_inds = 'results/local_ancestry/{model_name}/{sim_name}/{anal_name}.target_inds.txt',
	params:
		base_path = 'results/local_ancestry/{model_name}/{sim_name}/{anal_name}',
		chrom_id = lambda w: units.loc[(w.sim_name, w.anal_name)].chr,
		nind_ref = lambda w: units.loc[(w.sim_name, w.anal_name)].nind_ref,
		chr_len = lambda w: units.loc[w.sim_name, w.anal_name].chr_len,
	script:
		'../scripts/write_vcfs.py'


rule split_vcf:
	input:
		vcf =  'results/local_ancestry/{model_name}/{sim_name}/{anal_name}.genotypes.vcf.gz',
		samples = 'results/local_ancestry/{model_name}/{sim_name}/{anal_name}.{sample_group}.txt'
	output:
		vcf = 'results/local_ancestry/{model_name}/{sim_name}/{anal_name}.genotypes.{sample_group}.vcf.gz',
	params:
		bcftools = config['PATHS']['BCFTOOLS']
	shell:
		"""
		{params.bcftools} view --samples-file {input.samples} --output {output.vcf} --output-type z {input.vcf}
		{params.bcftools} index {output.vcf}
		"""

rule split_bcf:
	input:
		vcf =  'results/local_ancestry/{model_name}/{sim_name}/{anal_name}.genotypes.vcf.gz',
		samples = 'results/local_ancestry/{model_name}/{sim_name}/{anal_name}.{sample_group}.txt'
	output:
		bcf = 'results/local_ancestry/{model_name}/{sim_name}/{anal_name}.genotypes.{sample_group}.bcf',
	params:
		bcftools = config['PATHS']['BCFTOOLS']
	shell:
		"""
		{params.bcftools} view --samples-file {input.samples} --output {output.bcf} --output-type u {input.vcf}
		"""


rule true_local_ancestry:
	input:
		ancestry_ts = 'results/simulations/{model_name}/{sim_name}/{anal_name}.sample.tsz',
		site_ts = 'results/simulations/{model_name}/{sim_name}/{anal_name}.sample.filter.tsz',
		node_mapping = 'results/simulations/{model_name}/{sim_name}/{anal_name}.sample.filter.node_mapping.npy'
	output:
		tracts = 'results/local_ancestry/{model_name}/{sim_name}/{anal_name}.true_local_ancestry.tracts.hdf',
		site_matrix = 'results/local_ancestry/{model_name}/{sim_name}/{anal_name}.true_local_ancestry.site_matrix.npz',
		samples = 'results/local_ancestry/{model_name}/{sim_name}/{anal_name}.true_local_ancestry.samples.txt',
	params:
		admixture_time = lambda w: units.loc[(w.sim_name, w.anal_name)].admixture_time
	script:
		'../scripts/true_local_ancestry.py'


rule sample_sites:
	input:
		'results/simulations/{model_name}/{sim_name}/full.tsz'
	output:
		'results/simulations/{model_name}/{sim_name}/{anal_name}.sample.tsz'
	params:
		nind_admixed = lambda w: units.loc[(w.sim_name, w.anal_name)].nind_admixed,
		nind_ref = lambda w: units.loc[(w.sim_name, w.anal_name)].nind_ref,
		anal_seed = lambda w: units.loc[(w.sim_name, w.anal_name)].anal_seed,
		admixture_time = lambda w: units.loc[(w.sim_name, w.anal_name)].admixture_time
	script:
		"../scripts/sample.py"


rule filter_inds:
	input:
		'results/simulations/{model_name}/{sim_name}/{anal_name}.sample.tsz'
	output:
		tsz = 'results/simulations/{model_name}/{sim_name}/{anal_name}.sample.filter.tsz',
		node_mapping = 'results/simulations/{model_name}/{sim_name}/{anal_name}.sample.filter.node_mapping.npy'
	params:
		MAC_filter = lambda w: units.loc[(w.sim_name, w.anal_name)].MAC_filter,
		max_snps = lambda w: units.loc[(w.sim_name, w.anal_name)].max_snps,
		anal_seed = lambda w: units.loc[(w.sim_name, w.anal_name)].anal_seed
	script:
		"../scripts/filter.py"
