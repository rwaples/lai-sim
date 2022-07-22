"""Rules to evaluate the local ancestry methods."""
def all_sites(wildcards):
	import itertools
	sites = [f'results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/site.positions' for u in units.itertuples()],
	ret = list(itertools.chain(*sites))  # flatten
	return(ret)


def all_dosage(wildcards):
	import itertools
	flare = [f'results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/ancestry_dosage.flare.npz' for u in units.itertuples()],
	mosaic = [f'results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/ancestry_dosage.mosaic.npz' for u in units.itertuples()],
	rfmix2 = [f'results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/ancestry_dosage.rfmix2.npz' for u in units.itertuples()],
	ret = list(itertools.chain(*(flare + mosaic + rfmix2)))  # flatten
	return(ret)


def all_R2(wildcards):
	import itertools
	flare = [f'results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/SUMMARY/R2_score.flare.ancestry.tsv' for u in units.itertuples()],
	mosaic = [f'results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/SUMMARY/R2_score.mosaic.ancestry.tsv' for u in units.itertuples()],
	rfmix2 = [f'results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/SUMMARY/R2_score.rfmix2.ancestry.tsv' for u in units.itertuples()],
	ret = list(itertools.chain(*(flare + mosaic + rfmix2)))  # flatten
	return(ret)


def all_Q(wildcards):
	import itertools
	flare = [f'results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/SUMMARY/RMSD.flare.tsv' for u in units.itertuples()],
	mosaic = [f'results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/SUMMARY/RMSD.mosaic.tsv' for u in units.itertuples()],
	rfmix2 = [f'results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/SUMMARY/RMSD.rfmix2.tsv' for u in units.itertuples()],
	ret = list(itertools.chain(*(flare + mosaic + rfmix2)))  # flatten
	return(ret)


def all_runtime_benchmark(wildcards):
	import itertools
	flare = [f'results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/benchmark/run_flare.tsv' for u in units.itertuples()],
	mosaic = [f'results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/benchmark/run_mosaic.tsv' for u in units.itertuples()],
	rfmix2 = [f'results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/benchmark/run_RFMix2.tsv' for u in units.itertuples()],
	ret = list(itertools.chain(*(flare + mosaic + rfmix2)))  # flatten
	return(ret)


rule get_nsites:
	input:
		all_sites
	output:
		'results/reports/sites_report.txt',
	script:
		'../scripts/write_sites_report.py'


rule get_runtime_benchmark:
	input:
		all_runtime_benchmark
	output:
		'results/reports/runtime_report.txt',
	script:
		'../scripts/write_runtime_report.py'


rule make_R2_report:
	input:
		all_R2
	output:
		'results/reports/R2_report.txt',
	script:
		'../scripts/write_R2_report.py'


rule make_Q_report:
	input:
		all_Q
	output:
		'results/reports/Q_report.txt',
	script:
		'../scripts/write_Q_report.py'


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
	# conda:
	# 	'./env'
	script:
		"../scripts/get_Q_score.py"


rule write_qq_reports:
	input:
		true_la = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/ancestry_dosage.true.npz',
		flare_la = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/ancestry_dosage.flare.npz',
		mosaic_la = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/ancestry_dosage.mosaic.npz',
		rfmix2_la = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/ancestry_dosage.rfmix2.npz',

	output:
		flare_report = report('results/{model_name}/{sim_name}/{asc_name}/{anal_name}/DIAGNOSTICS/qq_flare.txt'),
		mosaic_report = report('results/{model_name}/{sim_name}/{asc_name}/{anal_name}/DIAGNOSTICS/qq_mosaic.txt'),
		rfmix_report = report('results/{model_name}/{sim_name}/{asc_name}/{anal_name}/DIAGNOSTICS/qq_rfmix.txt'),
	params:
		bcftools = config['PATHS']['BCFTOOLS'],
		nsource = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].nsource,
	script:
		"../scripts/write_qq_reports.py"


rule combine_qq_reports:
	input:
		reports = [f"results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/DIAGNOSTICS/qq_{{method}}.txt"
													for u in units.itertuples()],
	output:
		"results/reports/QQ.{method}.txt"
	params:
		reports = [f"results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/DIAGNOSTICS/qq_{{method}}.txt"
													for u in units.itertuples()],
	script:
		'../scripts/combine_qq_reports.py'
