rule reports:
	input:
		R2 = "results/reports/R2_report.txt",
		Q = "results/reports/Q_report.txt",
		QQ_bmix = "results/reports/QQ.bmix.txt",
		QQ_rfmix = "results/reports/QQ.rfmix.txt",
		QQ_mosaic = "results/reports/QQ.mosaic.txt",


rule ancestry_dosage_plots:
	input:
		bmix = [f"results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/DIAGNOSTICS/ancestry_dosage.{i}.bmix.pdf"
										for u in units.itertuples() for i in range(3)],
		rfmix2 = [f"results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/DIAGNOSTICS/ancestry_dosage.{i}.rfmix2.pdf"
												for u in units.itertuples() for i in range(3)],
		mosaic = [f"results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/DIAGNOSTICS/ancestry_dosage.{i}.mosaic.pdf"
												for u in units.itertuples() for i in range(3)]


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


rule plot_ancestry_dosage:
	input:
		true_la = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/true_local_ancestry.site_matrix.npz',
		mosaic_la = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/MOSAIC/la_probs.npz',
		rfmix2_la = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/RFMix2/rfmix2.fb.tsv.gz',
		bmix_la = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/bmix/bmix.anc.vcf.gz',
		sites_file = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/site.positions',
	output:
		bmix = report(
			expand(
				"results/{model_name}/{sim_name}/{asc_name}/{anal_name}/DIAGNOSTICS/ancestry_dosage.{i}.bmix.pdf",
				i=range(3),
				allow_missing=True
			)
		),
		rfmix2 = report(
			expand(
				"results/{model_name}/{sim_name}/{asc_name}/{anal_name}/DIAGNOSTICS/ancestry_dosage.{i}.rfmix2.pdf",
				i=range(3),
				allow_missing=True
			)
		),
		mosaic = report(
			expand(
				"results/{model_name}/{sim_name}/{asc_name}/{anal_name}/DIAGNOSTICS/ancestry_dosage.{i}.mosaic.pdf",
				i=range(3),
				allow_missing=True
			)
		),
	params:
		bcftools = config['PATHS']['BCFTOOLS'],
		nsource = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].nsource,
		path = "results/{model_name}/{sim_name}/{asc_name}/{anal_name}/DIAGNOSTICS/ancestry_dosage",
		format = 'pdf'
	script:
		'../scripts/plot_ancestry_dosage.py'
