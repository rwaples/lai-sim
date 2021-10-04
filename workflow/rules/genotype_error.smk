configfile: "./config.yml"

rule add_geno_error:
	input:
		gt = 'results/local_ancestry/{model_name}/{sim_name}/{anal_name}.genotypes.vcf.gz',
	output:
		vcf = 'results/local_ancestry/{model_name}/{sim_name}/{anal_name}.genotypes.witherror.vcf.gz',
		index = 'results/local_ancestry/{model_name}/{sim_name}/{anal_name}.genotypes.witherror.vcf.gz.csi',
	params:
		bcftools = config['PATHS']['BCFTOOLS'],
		gt_err = lambda w: units.loc[(w.sim_name, w.anal_name)].gt_err,
	shell:
		"""
		java -jar /projects/browning/software/add-err.jar {input.gt_err} | gzip > {output.vcf}
		{params.bcftools} index {output.vcf}
		"""
