"""Rules to fill in missing results from local ancestry analyses."""

rule fill_missing:
	input:
		reports = [f"results/{u.model_name}/{u.sim_name}/{u.asc_name}/{u.anal_name}/missing.LA.txt"
																for u in units.itertuples()],

rule add_missing_dosage:
	output:
		'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/missing.LA.txt',
	params:
		bmix = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/ancestry_dosage.bmix.npz',
		mosaic = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/ancestry_dosage.mosaic.npz',
		rfmix2 = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/ancestry_dosage.rfmix2.npz',
		log_bmix = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/benchmark/run_bmix.tsv',
		log_mosaic = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/benchmark/run_mosaic.tsv',
		log_rfmix2 = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/benchmark/run_RFMix2.tsv',

	log:
		'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/missing.LA.log',
	shell:
		"""
		touch {output}

		touch {params.bmix} {params.mosaic} {params.rfmix2} {params.log_bmix} {params.log_mosaic} {params.log_rfmix2}
		"""
