rule get_R2_score:
	input:
		true_dosage = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/ancestry_dosage.true.npz',
		target_dosage = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/ancestry_dosage.{method}.npz',
	output:
		R2_anc = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/SUMMARY/R2_score.{method}.ancestry.tsv',
		R2_ind = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/SUMMARY/R2_score.{method}.individuals.tsv',
	params:
		nsource = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].nsource,
	log:
		'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/SUMMARY/get_R2_score.{method}.log',
	script:
		"../scripts/get_R2_score.py"


def all_dosage(wildcards):
	import itertools
	bmix = [f'results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/ancestry_dosage.bmix.npz' for u in units.itertuples()],
	mosaic = [f'results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/ancestry_dosage.mosaic.npz' for u in units.itertuples()],
	rfmix2 = [f'results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/ancestry_dosage.rfmix2.npz' for u in units.itertuples()],
	ret = list(itertools.chain(*(bmix + mosaic + rfmix2)))  # flatten
	return(ret)


def all_R2(wildcards):
	import itertools
	bmix = [f'results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/SUMMARY/R2_score.bmix.ancestry.tsv' for u in units.itertuples()],
	mosaic = [f'results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/SUMMARY/R2_score.mosaic.ancestry.tsv' for u in units.itertuples()],
	rfmix2 = [f'results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/SUMMARY/R2_score.rfmix2.ancestry.tsv' for u in units.itertuples()],
	ret = list(itertools.chain(*(bmix + mosaic + rfmix2)))  # flatten
	return(ret)


def all_Q(wildcards):
	import itertools
	bmix = [f'results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/SUMMARY/RMSD.bmix.tsv' for u in units.itertuples()],
	mosaic = [f'results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/SUMMARY/RMSD.mosaic.tsv' for u in units.itertuples()],
	rfmix2 = [f'results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/SUMMARY/RMSD.rfmix2.tsv' for u in units.itertuples()],
	ret = list(itertools.chain(*(bmix + mosaic + rfmix2)))  # flatten
	return(ret)


rule make_R2_report_local:
	input:
		all_R2
	output:
		'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/SUMMARY/R2_report.txt',
	script:
		'../scripts/write_local_report.py'


rule make_Q_report_local:
	input:
		all_Q
	output:
		'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/SUMMARY/Q_report.txt',
	script:
		'../scripts/write_local_report.py'



rule make_R2_report_global:
	input:
		[f'results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/SUMMARY/R2_report.txt' for u in units.itertuples()]
	output:
		"results/reports/R2_report.txt",
	script:
		'../scripts/make_R2_report.py'


rule get_Q_score:
	input:
		true_dosage = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/ancestry_dosage.true.npz',
		target_path = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/ancestry_dosage.{method}.npz',
	output:
		Q_path = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/SUMMARY/Q.{method}.tsv',
		RMSD_path = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/SUMMARY/RMSD.{method}.tsv',
	params:
		nsource = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].nsource,
		Q_true = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/SUMMARY/Q.true.tsv',
	log:
		'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/SUMMARY/get_Q_score.{method}.log',
	conda:
		'./env'
	script:
		"../scripts/get_Q_score.py"


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