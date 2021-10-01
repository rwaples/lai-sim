configfile: "./config.yml"


rule add_geno_error:
	input:
		[f'results/local_ancestry/{u.model_name}/{u.sim_name}/{u.anal_name}.genotypes.witherror.vcf.gz' for u in units.itertuples()]

rule run_add_err:
	input:
		gt = 'results/local_ancestry/{model_name}/{sim_name}/{anal_name}.genotypes.vcf.gz',
	output:
		vcf = 'results/local_ancestry/{model_name}/{sim_name}/{anal_name}.genotypes.witherror.vcf.gz',
		index = 'results/local_ancestry/{model_name}/{sim_name}/{anal_name}.genotypes.witherror.vcf.gz.csi',
	params:
		prefix = 'results/local_ancestry/{model_name}/{sim_name}/{anal_name}.phased.target_inds',
		bcftools = config['PATHS']['BCFTOOLS'],
		seed = lambda w: units.loc[(w.sim_name, w.anal_name)].anal_seed,
		gt_err = lambda w: units.loc[(w.sim_name, w.anal_name)].gt_err,
	shell:
		"""
		java -jar /projects/browning/software/add-err.jar {input.gt_err} | gzip > {output.vcf}
		{params.bcftools} index {output.vcf}
		"""
