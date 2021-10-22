
# change this to use simpler 3-wild names
rule export_mosaic:
	input:
		la_results = 'results/{model_name}/{sim_name}/{anal_name}/MOSAIC/localanc_admixed_{nsource}way_1-{naming_mess}.RData',
		model_results = 'results/{model_name}/{sim_name}/{anal_name}/MOSAIC/admixed_{nsource}way_1-{naming_mess}.RData',
	output:
		path = 'results/{model_name}/{sim_name}/{anal_name}/MOSAIC/la_probabilites.{nsource}way_1-{naming_mess}.RData',
	params:
		input_dir = 'results/{model_name}/{sim_name}/{anal_name}/MOSAIC/input/',
		simple_output = 'results/{model_name}/{sim_name}/{anal_name}/MOSAIC/la_rename.RData',
	script:
		"../scripts/export_mosaic_results.R"





# these functions are meant to predict the naming scheme used by MOSAIC
# but these sometime fails to predict the correct name. This seems to be due
# to some aspects being determined at runtime, such as max.donors
# see workflow/rules/common.smk for the predctions of nsource (easy) and (naming_mess) hard
def localanc_name(wildcards):
	#model_name = wildcards.model_name
	#sim_name = wildcards.sim_name
	#anal_name = wildcards.anal_name
	#nsource = units.loc[(sim_name, anal_name)].nsource
	#naming_mess = units.loc[(sim_name, anal_name)].naming_mess
	#return f'results/{model_name}/{sim_name}/{anal_name}/MOSAIC/localanc_admixed_{nsource}way_1-{naming_mess}.RData'
	return []

def admixed_name(wildcards):
	#model_name = wildcards.model_name
	#sim_name = wildcards.sim_name
	#anal_name = wildcards.anal_name
	#nsource = units.loc[(sim_name, anal_name)].nsource
	#naming_mess = units.loc[(sim_name, anal_name)].naming_mess
	#return f'results/{model_name}/{sim_name}/{anal_name}/MOSAIC/admixed_{nsource}way_1-{naming_mess}.RData'
	return []

def aggregate_input(wildcards):
	'''
	aggregate the file names of the random number of files
	this is now just serving to defer
	'''
	#checkpoint_output = checkpoints.run_mosaic.get(**wildcards).output[0]
	return []



rule export_mosaic:
	input:
		la_results = 'results/{model_name}/{sim_name}/{anal_name}/MOSAIC/localanc_admixed_{nsource}way_1-{naming_mess}.RData',
		model_results = 'results/{model_name}/{sim_name}/{anal_name}/MOSAIC/admixed_{nsource}way_1-{naming_mess}.RData',
	output:
		path = 'results/{model_name}/{sim_name}/{anal_name}/MOSAIC/la_probabilites.{nsource}way_1-{naming_mess}.RData',
	params:
		input_dir = 'results/{model_name}/{sim_name}/{anal_name}/MOSAIC/input/',
		simple_output = 'results/{model_name}/{sim_name}/{anal_name}/MOSAIC/la_rename.RData',
	script:
		"../scripts/export_mosaic_results.R"



checkpoint dummy:
	input:
	output: touch('results/{model_name}/{sim_name}/{anal_name}/MOSAIC/.dummy')



def generate_mosaic_name(wildcards):
	import glob
	print(wildcards)
	check = checkpoints.export_mosaic.get(**wildcards)
	model_name = wildcards.model_name
	sim_name = wildcards.sim_name
	anal_name = wildcards.anal_name
	nsource = check.nsource
	naming_mess = check.naming_mess
	paths = glob.glob(f'results/{model_name}/{sim_name}/{anal_name}/MOSAIC/la_probabilites.*.RData')
	print(paths)
	assert len(paths)==1
	return paths[0]


rule mosiac_dummy:
	input:
		generate_mosaic_name
	output:
		'results/{model_name}/{sim_name}/{anal_name}/MOSAIC/la_probs.RData',
	shell:
		"cp {input} {output}"


rule run_mosaic:
	input:
		'programs/MOSAIC/MOSAIC/mosaic.R',
		'results/{model_name}/{sim_name}/{anal_name}/MOSAIC/input/admixedgenofile.22',
	output:
		#la_results = "results/{model_name}/{sim_name}/{anal_name}/MOSAIC/localanc_admixed_{nsource}way_1-{naming_mess}.RData",
		#model_results = "results/{model_name}/{sim_name}/{anal_name}/MOSAIC/admixed_{nsource}way_1-{naming_mess}.RData",
		la_results = 'results/{model_name}/{sim_name}/{anal_name}/MOSAIC/localanc_admixed.RData',
		model_results = 'results/{model_name}/{sim_name}/{anal_name}/MOSAIC/admixed.RData',
	log:
		#'results/{model_name}/{sim_name}/{anal_name}/MOSAIC/mosaic.{nsource}way_1-{naming_mess}.log'
		'results/{model_name}/{sim_name}/{anal_name}/MOSAIC/mosaic.log',
	benchmark:
		#'results/{model_name}/{sim_name}/{anal_name}/benchmark/run_mosaic.{nsource}way_1-{naming_mess}.tsv',
		'results/{model_name}/{sim_name}/{anal_name}/benchmark/run_mosaic.tsv',
	params:
		mosaic = config['PATHS']['MOSAIC'],
		input_folder = 'results/{model_name}/{sim_name}/{anal_name}/MOSAIC/input/',
		base_folder = 'results/{model_name}/{sim_name}/{anal_name}/MOSAIC',
		seed = lambda w: units.loc[(w.sim_name, w.anal_name)].anal_seed,
		nsource = lambda w: units.loc[(w.sim_name, w.anal_name)].nsource,
		nthreads = lambda w: units.loc[(w.sim_name, w.anal_name)].nthreads,
		#rename_la = 'results/{model_name}/{sim_name}/{anal_name}/MOSAIC/localanc_admixed.RData',
		#rename_model = 'results/{model_name}/{sim_name}/{anal_name}/MOSAIC/admixed.RData',
	shell:
		"""

		Rscript {params.mosaic} --seed {params.seed} --maxcores {params.nthreads} --chromosomes 22:22 --ancestries {params.nsource} admixed {params.input_folder} 2>&1 | tee {log}

		"""

		"""

		# move the results into the MOSAIC folder
		mv MOSAIC_RESULTS/* {params.base_folder}
		mv MOSAIC_PLOTS/* {params.base_folder}
		mv FREQS/* {params.base_folder}
		rmdir MOSAIC_RESULTS
		rmdir MOSAIC_PLOTS
		rmdir FREQS

		# try to fix the naming mess here
		cd {params.base_folder}
		find localanc_admixed_*.RData -exec cp {} localanc_admixed.RData \;
		find admixed_*.RData -exec cp {} admixed.RData \;

		"""
