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


rule write_qq_reports:
	input:
		true_la = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/true_local_ancestry.site_matrix.npz',
		mosaic_la = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/MOSAIC/la_probs.npz',
		rfmix2_la = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/RFMix2/rfmix2.fb.tsv.gz',
		bmix_la = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/bmix/bmix.anc.vcf.gz',
		sites_file = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/site.positions',
	output:
		bmix_report = report('results/{model_name}/{sim_name}/{asc_name}/{anal_name}/DIAGNOSTICS/qq_bmix.txt'),
		mosaic_report = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/DIAGNOSTICS/qq_mosaic.txt',
		rfmix_report = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/DIAGNOSTICS/qq_rfmix.txt',
	params:
		bcftools = config['PATHS']['BCFTOOLS'],
		nsource = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].nsource,
	script:
		"../scripts/write_qq_reports.py"


rule plot_qq_reports:
	input:
		bmix_report = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/DIAGNOSTICS/qq_bmix.txt',
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
