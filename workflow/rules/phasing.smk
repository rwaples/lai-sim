configfile: "./config.yml"


rule phasing:
	input:
		[f'results/local_ancestry/{u.model_name}/{u.sim_name}/{u.anal_name}.phased.target_inds.vcf.gz' for u in units.itertuples()],
		[f'results/local_ancestry/{u.model_name}/{u.sim_name}/{u.anal_name}.phased.reference_inds.vcf.gz' for u in units.itertuples()]

rule run_beagle:
	input:
		gt = 'results/local_ancestry/{model_name}/{sim_name}/{anal_name}.genotypes.witherror.vcf.gz',
		map = 'results/local_ancestry/{model_name}/{sim_name}/{anal_name}.plink_map.txt',
	output:
		vcf = 'results/local_ancestry/{model_name}/{sim_name}/{anal_name}.phased.vcf.gz',
		index = 'results/local_ancestry/{model_name}/{sim_name}/{anal_name}.phased.vcf.gz.csi',
	log:
		'results/local_ancestry/{model_name}/{sim_name}/{anal_name}.phased.vcf.beagle.log',
	params:
		prefix = 'results/local_ancestry/{model_name}/{sim_name}/{anal_name}.phased',
		bcftools = config['PATHS']['BCFTOOLS'],
		seed = lambda w: units.loc[(w.sim_name, w.anal_name)].anal_seed,
	shell:
		"java -jar programs/BEAGLE/beagle.29May21.d6d.jar "
		"gt={input.gt} "
		"map={input.map} "
		"seed={params.seed} "
		"out={params.prefix} 2>&1 | tee {log} "

		"""

		{params.bcftools} index {output.vcf}
		"""
