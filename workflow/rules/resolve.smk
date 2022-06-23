"""Rules to fill in missing results from local ancestry analyses."""
rule fill_missing:
	input:
		reports = [f"results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/missing.LA.txt"
													for u in units.itertuples()],

rule add_missing_dosage:
	output:
		'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/missing.LA.txt',
	params:
		flare = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/ancestry_dosage.flare.npz',
		mosaic = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/ancestry_dosage.mosaic.npz',
		rfmix2 = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/ancestry_dosage.rfmix2.npz',
		log_flare = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/benchmark/run_flare.tsv',
		log_mosaic = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/benchmark/run_mosaic.tsv',
		log_rfmix2 = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/benchmark/run_RFMix2.tsv',

	log:
		'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/missing.LA.log',
	shell:
		"""
		touch {output}

		touch {params.flare} {params.mosaic} {params.rfmix2} {params.log_flare} {params.log_mosaic} {params.log_rfmix2}
		"""
