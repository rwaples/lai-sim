rule ancestry_dosage_plots:
	input:
		flare = [f"results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/DIAGNOSTICS/ancestry_dosage.{i}.flare.pdf"
										for u in units.itertuples() for i in range(3)],
		rfmix2 = [f"results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/DIAGNOSTICS/ancestry_dosage.{i}.rfmix2.pdf"
												for u in units.itertuples() for i in range(3)],
		mosaic = [f"results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/DIAGNOSTICS/ancestry_dosage.{i}.mosaic.pdf"
												for u in units.itertuples() for i in range(3)]


rule plot_ancestry_dosage:
	input:
		true_la = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/true_local_ancestry.site_matrix.npz',
		mosaic_la = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/MOSAIC/la_probs.npz',
		rfmix2_la = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/RFMix2/rfmix2.fb.tsv.gz',
		flare_la = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/flare/flare.anc.vcf.gz',
		sites_file = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/site.positions',
	output:
		flare = report(
			expand(
				"results/{model_name}/{sim_name}/{asc_name}/{anal_name}/DIAGNOSTICS/ancestry_dosage.{i}.flare.pdf",
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


rule plot_pop_coal_time:
	input:
		site_ts = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/sample.filter.tsz'
	output:
		plot = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/DIAGNOSTICS/pop_coal_time.png'
	script:
		'../scripts/plot_pop_coal_time.py'


rule plot_pairwise_Fst:
	input:
		site_ts = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/sample.filter.tsz',
	output:
		plot = report(
			'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/DIAGNOSTICS/pairwise_Fst.png',
			caption="../report/pairwise_Fst.rst",
			category="Diagnostics",
			subcategory="{model_name}  {sim_name}  {asc_name}  {anal_name}"
		),
	script:
		'../scripts/plot_pairwise_Fst.py'


rule plot_local_ancestry:
	input:
		la_mat = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/true_local_ancestry.site_matrix.npz',
	output:
		plot = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/DIAGNOSTICS/true_local_ancestry.png',
	script:
		'../scripts/plot_true_local_ancestry.py'


rule plot_qq_reports:
	input:
		flare_report = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/DIAGNOSTICS/qq_flare.txt',
		mosaic_report = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/DIAGNOSTICS/qq_mosaic.txt',
		rfmix_report = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/DIAGNOSTICS/qq_rfmix.txt',
	output:
		plot_path = report(
			'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/DIAGNOSTICS/qq.png',
			caption="../report/qq_plot.rst",
			category="Diagnostics",
			subcategory="{model_name}  {sim_name}  {asc_name}  {anal_name}"
		),
	script:
		'../scripts/plot_qq.py'
