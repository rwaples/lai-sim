"""Create data and files for LA analyses, prior to phasing"""
rule write_sites:
	input:
		bcf = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/genotypes.bcf',
	output:
		sites = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/site.positions',
	params:
		BCFTOOLS = config['PATHS']['BCFTOOLS'],
	shell:
		'{params.BCFTOOLS} query -f "%POS\\n" {input.bcf} > {output.sites}'


rule write_vcfs:
	input:
		site_matrix = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/true_local_ancestry.site_matrix.npz',
		site_ts = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/sample.filter.tsz',
	output:
		phased_ancestry_vcf = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/la_true.vcf.gz',
		phased_genotype_vcf = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/genotypes.vcf.gz',
		phased_genotype_bcf = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/genotypes.bcf',
		#genetic_map = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/genetic_map.txt',
		#plink_map = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/plink_map.txt',
		sample_map = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/sample_map.txt',
		reference_inds = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/reference_inds.txt',
		target_inds = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/target_inds.txt',
		#sites = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/site.positions',
	params:
		base_path = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}',
		chrom_id = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].chr,
		nind_ref = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].nind_ref,
		# chr_len = lambda w: units.loc[w.sim_name, w.asc_name, w.anal_name].chr_len,
		target_pop = lambda w: units.loc[w.sim_name, w.asc_name, w.anal_name].target_pop,
	script:
		'../scripts/write_vcfs.py'


rule split_vcf:
	input:
		vcf = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/phased.vcf.gz',
		samples = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/{sample_group}.txt',
	output:
		vcf = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/phased.{sample_group}.vcf.gz',
	params:
		BCFTOOLS = config['PATHS']['BCFTOOLS']
	shell:
		"""
		{params.BCFTOOLS} view --samples-file {input.samples} --output-file {output.vcf} --output-type z {input.vcf}
		{params.BCFTOOLS} index {output.vcf}
		"""


rule split_bcf:
	input:
		vcf = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/phased.vcf.gz',
		samples = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/{sample_group}.txt'
	output:
		bcf = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/phased.{sample_group}.bcf',
	params:
		bcftools = config['PATHS']['BCFTOOLS']
	shell:
		"""
		{params.bcftools} view --samples-file {input.samples} --output-file {output.bcf} --output-type b {input.vcf}
		{params.bcftools} index {output.bcf}
		"""


rule true_local_ancestry:
	input:
		ancestry_ts = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/sample.tsz',
		site_ts = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/sample.filter.tsz',
		node_mapping = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/sample.filter.node_mapping.npy'
	output:
		tracts = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/true_local_ancestry.tracts.hdf',
		site_matrix = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/true_local_ancestry.site_matrix.npz',
		samples = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/true_local_ancestry.samples.txt',
	params:
		admixture_time = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].admixture_time,
		target_pop = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].target_pop,
		nind_admixed = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].nind_admixed,
	script:
		'../scripts/true_local_ancestry.py'


rule ascertain:
	input:
		'results/{model_name}/{sim_name}/full.tsz'
	output:
		'results/{model_name}/{sim_name}/{asc_name}/ascertained.tsz'
	log:
		'results/{model_name}/{sim_name}/{asc_name}/ascertained.log'
	params:
		asc_MAC = lambda w: units.loc[(w.sim_name, w.asc_name)].asc_MAC[0],
		asc_MAF = lambda w: units.loc[(w.sim_name, w.asc_name)].asc_MAF[0],
		asc_maxsites = lambda w: units.loc[(w.sim_name, w.asc_name)].asc_maxsites[0],
		asc_seed = lambda w: units.loc[(w.sim_name, w.asc_name)].asc_seed[0],
		target_pop = lambda w: units.loc[(w.sim_name, w.asc_name)].target_pop[0],
	script:
		"../scripts/ascertain.py"


rule sample_sites:
	input:
		'results/{model_name}/{sim_name}/{asc_name}/ascertained.tsz'
	output:
		'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/sample.tsz'
	params:
		nind_admixed = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].nind_admixed,
		nind_ref = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].nind_ref,
		anal_seed = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].anal_seed,
		admixture_time = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].admixture_time,
	script:
		"../scripts/sample.py"


rule filter_inds:
	input:
		'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/sample.tsz'
	output:
		tsz = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/sample.filter.tsz',
		node_mapping = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/sample.filter.node_mapping.npy'
	params:
		MAC_filter = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].MAC_filter,
		max_snps = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].max_snps,
		anal_seed = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].anal_seed,
	script:
		"../scripts/filter_sites.py"


rule add_geno_error:
	input:
		gt = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/genotypes.vcf.gz',
	output:
		vcf = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/genotypes.witherror.vcf.gz',
		index = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/genotypes.witherror.vcf.gz.csi',
	params:
		bcftools = config['PATHS']['BCFTOOLS'],
		ADD_ERR = config['PATHS']['ADD_ERR'],
		gt_err = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].gt_err,
	shell:
		"""
		zcat {input.gt} | java -jar {params.ADD_ERR} {params.gt_err} | bgzip > {output.vcf}
		{params.bcftools} index {output.vcf}
		"""
