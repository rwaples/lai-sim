rule sims_done:
	input:
		[f'results/{s.model_name}/{s.sim_name}/full.tsz' for s in simulations.itertuples()]


rule recap_and_mutate:
	input:
		"results/{model_name}/{sim_name}/from_slim.trees",
	output:
		'results/{model_name}/{sim_name}/full.tsz'
	params:
		# params are used by the script
		ancestral_Ne = lambda w: simulations.loc[w.sim_name].ancestral_Ne,
		chr = lambda w: simulations.loc[w.sim_name].chr,
		chr_len = lambda w: simulations.loc[w.sim_name].chr_len,
		mutation_rate = lambda w: simulations.loc[w.sim_name].mutation_rate,
		sim_seed = lambda w: simulations.loc[w.sim_name].sim_seed,
	script:
		"../scripts/recap_and_mutate.py"


rule simulate_admixture:
	output:
		"results/{model_name}/{sim_name}/from_slim.trees"
	log:
		"results/{model_name}/{sim_name}/from_slim.trees.log"
	params:
		slim_path = config["PATHS"]["SLiM"],
		sim_seed = lambda w: simulations.loc[w.sim_name].sim_seed,
		slim_script = lambda w: simulations.loc[w.sim_name].slim_script_path,
	shell:
		"""{params.slim_path} -seed {params.sim_seed} -d 'trees_file="{output}"' {params.slim_script} 2>&1 | tee {log}"""


rule slim_scripts_present:
	input:
		[s.slim_script_path for s in simulations.itertuples()]
