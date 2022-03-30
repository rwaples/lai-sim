
rule get_Q_score:
	input:
		true_la = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/true_local_ancestry.site_matrix.npz',
		mosaic_la = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/MOSAIC/la_probs.npz',
		rfmix2_la = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/RFMix2/rfmix2.fb.tsv.gz',
		bmix_la = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/bmix/bmix.anc.vcf.gz',
		sites_file = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/site.positions',
	output:
		RMSD_path = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/SUMMARY/Q_RMSD.tsv',
		Q_true_path = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/SUMMARY/Q_true.tsv',
		Q_bmix_path = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/SUMMARY/Q_bmix.tsv',
		Q_mosaic_path = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/SUMMARY/Q_mosaic.tsv',
		Q_rfmix_path = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/SUMMARY/Q_rfmix.tsv',
	params:
		bcftools = config['PATHS']['BCFTOOLS'],
		nsource = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].nsource,
	log:
		'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/SUMMARY/get_Q_score.log',
	conda:
		'./env'
	script:
		"../scripts/get_Q_score.py"


rule get_R2_score:
	input:
		true_la = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/true_local_ancestry.site_matrix.npz',
		mosaic_la = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/MOSAIC/la_probs.npz',
		rfmix2_la = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/RFMix2/rfmix2.fb.tsv.gz',
		bmix_la = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/bmix/bmix.anc.vcf.gz',
		sites_file = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/site.positions',
	output:
		R2_anc = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/SUMMARY/R2_score.ancestry.tsv',
		R2_ind = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/SUMMARY/R2_score.individuals.tsv',
	params:
		bcftools = config['PATHS']['BCFTOOLS'],
		nsource = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].nsource,
	log:
		'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/SUMMARY/get_R2_score.log',
	script:
		"../scripts/get_R2_score.py"


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
