def gather_R2(wildcards):
	checkpoints.get_R2_score.get(**wildcards)
	#methods = glob_wildcards(f'results/{wildcards.model_name}/{wildcards.sim_name}/{wildcards.asc_name}/{wildcards.anal_name}/SUMMARY/R2_score.{{method}}.ancestry.tsv').method
	#ret = expand('results/{wildcards.model_name}/{wildcards.sim_name}/{wildcards.asc_name}/{wildcards.anal_name}/SUMMARY/R2_score.{{method}}.ancestry.tsv', method=methods)
	#meth, = glob_wildcards(
	#	f'results/{wildcards.model_name}/{wildcards.sim_name}/{wildcards.asc_name}/{wildcards.anal_name}/SUMMARY/R2_score.{{method}}.ancestry.tsv'
	#).method
	ret = expand(
		'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/SUMMARY/R2_score.{method}.ancestry.tsv',
		model_name=wildcards.model_name,
		sim_name=wildcards.sim_name,
		asc_name=wildcards.asc_name,
		anal_name=wildcards.anal_name,
		method=method
	)
	return(ret)


#rule make_R2_report_global:
#	input:
#		gather_R2
#	output:
#		"results/reports/R2_report.txt",
#	params:
#		dummy = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/ancestry_dosage.true.npz',
#		# files = [f"results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/SUMMARY/R2_score.ancestry.tsv"
#		# 									for u in units.itertuples()],
#		# names = [f"{u.model_name}\t{u.sim_name}\t{u.asc_name}\t{u.anal_name}" for u in units.itertuples()],
#		# nind_ref = [units.loc[(u.sim_name, u.asc_name, u.anal_name)].nind_ref for u in units.itertuples()],
#	script:
#		'../scripts/make_R2_report.py'
