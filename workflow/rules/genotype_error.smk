configfile: "./config.yml"

rule add_geno_error:
	input:
		gt = 'results/{model_name}/{sim_name}/{anal_name}/genotypes.vcf.gz',
	output:
		vcf = 'results/{model_name}/{sim_name}/{anal_name}/genotypes.witherror.vcf.gz',
		index = 'results/{model_name}/{sim_name}/{anal_name}/genotypes.witherror.vcf.gz.csi',
	params:
		bcftools = config['PATHS']['BCFTOOLS'],
		ADD_ERR = config['PATHS']['ADD_ERR'],
		gt_err = lambda w: units.loc[(w.sim_name, w.anal_name)].gt_err,
	shell:
		"""
		zcat {input.gt} | java -jar {params.ADD_ERR} {params.gt_err} | bgzip > {output.vcf}
		{params.bcftools} index {output.vcf}
		"""
