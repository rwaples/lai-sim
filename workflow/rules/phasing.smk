rule run_beagle:
	input:
		gt = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/genotypes.witherror.vcf.gz',
		map = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/plink_map.txt',
	output:
		vcf = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/phased.vcf.gz',
		index = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/phased.vcf.gz.csi',
	log:
		'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/phased.vcf.beagle.log',
	params:
		prefix = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/phased',
		bcftools = config['PATHS']['BCFTOOLS'],
		BEAGLE = config['PATHS']['BEAGLE'],
		seed = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].anal_seed,
	shell:
		"java -jar {params.BEAGLE} "
		"gt={input.gt} "
		"map={input.map} "
		"seed={params.seed} "
		"out={params.prefix} 2>&1 | tee {log} "

		"""

		{params.bcftools} index {output.vcf}
		"""
