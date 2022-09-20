wildcard_constraints:
	anal_name = '\w+',
	asc_name = '\w+',
	sim_name = '\w+',


rule simulate_admixture:
	input:
		'results/{model_name}/{sim_name}/slim_map.txt',
	output:
		temp('results/{model_name}/{sim_name}/from_slim.trees'),
	log:
		'results/{model_name}/{sim_name}/from_slim.trees.log'
	params:
		slim_path = config['PATHS']['SLiM'],
		sim_seed = lambda w: simulations.loc[w.sim_name].sim_seed,
		slim_script = lambda w: simulations.loc[w.sim_name].slim_script_path,
		max_bp = lambda w: simulations.loc[w.sim_name].max_bp,
		size_A = lambda w: simulations.loc[w.sim_name].size_A,
		size_B = lambda w: simulations.loc[w.sim_name].size_B,
		size_admixed = lambda w: simulations.loc[w.sim_name].size_admixed,
		time_split = lambda w: simulations.loc[w.sim_name].time_split,
	benchmark:
		'results/{model_name}/{sim_name}/benchmark/simulate_admixture.tsv',
	shell:
		"""{params.slim_path} \
			-seed {params.sim_seed} \
			-d 'slim_map="{input}"' \
			-d 'trees_file="{output}"' \
			-d 'max_bp="{params.max_bp}"' \
			-d 'sA="{params.size_A}"' \
			-d 'sB="{params.size_B}"' \
			-d 'sadmixed="{params.size_admixed}"' \
			-d 'tsplit="{params.time_split}"' \
			{params.slim_script} 2>&1 | tee {log}"""


rule recap_and_mutate:
	input:
		ts_path = 'results/{model_name}/{sim_name}/from_slim.trees',
	output:
		'results/{model_name}/{sim_name}/full.tsz'
	log:
		'results/{model_name}/{sim_name}/recap_and_mutate.log'
	params:
		ancestral_Ne = lambda w: simulations.loc[w.sim_name].ancestral_Ne,
		#chr = lambda w: simulations.loc[w.sim_name].chr,
		#chr_len = lambda w: simulations.loc[w.sim_name].chr_len,
		rec_map_path = lambda w: simulations.loc[w.sim_name].rec_map_path,
		max_bp = lambda w: simulations.loc[w.sim_name].max_bp,
		mutation_rate = lambda w: simulations.loc[w.sim_name].mutation_rate,
		sim_seed = lambda w: simulations.loc[w.sim_name].sim_seed,
		admixture_time = lambda w: simulations.loc[w.sim_name].admixture_time,
	benchmark:
		'results/{model_name}/{sim_name}/benchmark/recap_and_mutate.tsv',
	script:
		'../scripts/recap_and_mutate.py'


rule write_slim_map:
	output:
		slim_map_path = 'results/{model_name}/{sim_name}/slim_map.txt'
	params:
		rec_map_path = lambda w: simulations.loc[w.sim_name].rec_map_path,
		max_bp = lambda w: simulations.loc[w.sim_name].max_bp,
	script:
		'../scripts/write_slim_map.py'
