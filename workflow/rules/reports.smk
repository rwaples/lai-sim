rule report:
	input:
		[f"results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/SUMMARY/R2_score.ancestry.tsv"
			for u in units.itertuples()]
	output:
		"results/reports/R2_report.txt"
	params:
		files = [f"results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/SUMMARY/R2_score.ancestry.tsv"
			for u in units.itertuples()],
		names = [f"{u.model_name}\t{u.sim_name}\t{u.asc_name}\t{u.anal_name}" for u in units.itertuples()],
		nind_ref = [units.loc[(u.sim_name, u.asc_name, u.anal_name)].nind_ref for u in units.itertuples()],
	script:
		'../scripts/make_R2_report.py'
