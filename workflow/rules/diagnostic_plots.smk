rule plot_pop_coal_time:
	input:
		site_ts = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/sample.filter.tsz'
	output:
		plot = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/DIAGNOSTICS/pop_coal_time.png'
	script:
		'../scripts/plot_pop_coal_time.py'


rule plot_pairwise_Fst:
	input:
		site_ts = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/sample.filter.tsz'
	output:
		plot = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/DIAGNOSTICS/pairwise_Fst.png'
	script:
		'../scripts/plot_pairwise_Fst.py'


rule plot_local_ancestry:
	input:
		la_mat = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/true_local_ancestry.site_matrix.npz'
	output:
		plot = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/DIAGNOSTICS/true_local_ancestry.png'
	script:
		'../scripts/plot_true_local_ancestry.py'
