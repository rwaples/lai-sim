rule reports:
	input:
		R2 = "results/reports/R2_report.txt",
		Q = "results/reports/Q_report.txt",
		QQ_bmix = "results/reports/QQ.bmix.txt",
		QQ_rfmix = "results/reports/QQ.rfmix.txt",
		QQ_mosaic = "results/reports/QQ.mosaic.txt",


rule R2_report:
	input:
		[f"results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/SUMMARY/R2_score.ancestry.tsv"
			for u in units.itertuples()],
	output:
		"results/reports/R2_report.txt",
	params:
		files = [f"results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/SUMMARY/R2_score.ancestry.tsv"
											for u in units.itertuples()],
		names = [f"{u.model_name}\t{u.sim_name}\t{u.asc_name}\t{u.anal_name}" for u in units.itertuples()],
		nind_ref = [units.loc[(u.sim_name, u.asc_name, u.anal_name)].nind_ref for u in units.itertuples()],
	script:
		'../scripts/make_R2_report.py'


rule Q_report:
	input:
		[f"results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/SUMMARY/Q_RMSD.tsv"
			for u in units.itertuples()]
	output:
		"results/reports/Q_report.txt"
	params:
		files = [f"results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/SUMMARY/Q_RMSD.tsv"
											for u in units.itertuples()],
		names = [f"{u.model_name}\t{u.sim_name}\t{u.asc_name}\t{u.anal_name}" for u in units.itertuples()],
		nind_ref = [units.loc[(u.sim_name, u.asc_name, u.anal_name)].nind_ref for u in units.itertuples()],
	script:
		'../scripts/make_Q_report.py'


rule combine_qq_reports:
	input:
		reports = [f"results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/DIAGNOSTICS/qq_{{lai_program}}.txt"
													for u in units.itertuples()],
	output:
		"results/reports/QQ.{lai_program}.txt"
	params:
		reports = [f"results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/DIAGNOSTICS/qq_{{lai_program}}.txt"
													for u in units.itertuples()],
	script:
		'../scripts/combine_qq_reports.py'
