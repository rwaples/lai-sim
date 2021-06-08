configfile: "./config.yml"


rule phasing:
	input:
		[f'results/local_ancestry/{u.model_name}/{u.sim_name}/{u.anal_name}.phased.target_inds.vcf.gz' for u in units.itertuples()]

rule run_beagle:
	input:
		gt = 'results/local_ancestry/{model_name}/{sim_name}/{anal_name}.genotypes.target_inds.vcf.gz',
		ref = 'results/local_ancestry/{model_name}/{sim_name}/{anal_name}.genotypes.reference_inds.vcf.gz',
		map = 'results/local_ancestry/{model_name}/{sim_name}/{anal_name}.plink_map.txt',
	output:
		vcf = 'results/local_ancestry/{model_name}/{sim_name}/{anal_name}.phased.target_inds.vcf.gz',
		index = 'results/local_ancestry/{model_name}/{sim_name}/{anal_name}.phased.target_inds.vcf.gz.csi',
	log:
		'results/local_ancestry/{model_name}/{sim_name}/{anal_name}.phased.target_inds.vcf.beagle.log',
	params:
		prefix = 'results/local_ancestry/{model_name}/{sim_name}/{anal_name}.phased.target_inds',
		bcftools = config['PATHS']['BCFTOOLS']
	shell:
		"java -jar programs/BEAGLE/beagle.29May21.d6d.jar "
		"gt={input.gt} "
		"ref={input.ref} "
		"map={input.map} "
		"out={params.prefix} 2>&1 | tee {log}"

		"""

		{params.bcftools} index {output.vcf}
		"""
