configfile: "./config.yml"
wildcard_constraints:
	anal_name="\w+"


rule export_mosaic:
	input:
		la_results = 'results/{model_name}/{sim_name}/{anal_name}/MOSAIC/localanc_admixed_{nsource}way_1-{naming_mess}.RData',
		model_results = 'results/{model_name}/{sim_name}/{anal_name}/MOSAIC/admixed_{nsource}way_1-{naming_mess}.RData',
	output:
		path = 'results/{model_name}/{sim_name}/{anal_name}/MOSAIC/la_probabilites.{nsource}way_1-{naming_mess}.RData',
	params:
		input_dir = 'results/{model_name}/{sim_name}/{anal_name}/MOSAIC/input/',
	script:
		"../scripts/export_mosaic_results.R"


rule summarize_rfmix2:
	input:
		inferred_la = 'results/{model_name}/{sim_name}/{anal_name}/RFMix2/rfmix2.fb.tsv'
	output:
		diploid = 'results/{model_name}/{sim_name}/{anal_name}/RFMix2/diploid_la.hdf',
	params:
		threshold = 0.9
	script:
		"../scripts/summarize_rfmix2.py"


rule run_bmix:
	input:
		target_vcf = 'results/{model_name}/{sim_name}/{anal_name}/phased.target_inds.vcf.gz',
		reference_vcf = 'results/{model_name}/{sim_name}/{anal_name}/phased.reference_inds.vcf.gz',
		sample_map = 'results/{model_name}/{sim_name}/{anal_name}/sample_map.txt',
		genetic_map = 'results/{model_name}/{sim_name}/{anal_name}/plink_map.txt',
	output:
		'results/{model_name}/{sim_name}/{anal_name}/bmix/bmix.anc.vcf.gz',
	params:
		BMIX = config['PATHS']['BMIX'],
		bcftools = config['PATHS']['BCFTOOLS'],
		prefix = 'results/{model_name}/{sim_name}/{anal_name}/bmix/bmix',
		nthreads = 8,
		seed = lambda w: units.loc[(w.sim_name, w.anal_name)].anal_seed,
	shell:
		"java -Xmx30g -jar {params.BMIX} "
		"ref={input.reference_vcf} "
		"ref-panel={input.sample_map} "
		"gt={input.target_vcf} "
		"map={input.genetic_map} "
		"out={params.prefix} "
		"probs=true "
		"nthreads={params.nthreads} "
		"seed={params.seed} "
		"min-maf=0 "

		"""

		{params.bcftools} index {params.prefix}.vcf.gz
		"""




rule run_mosaic:
	input:
		'programs/MOSAIC/MOSAIC/mosaic.R',
		'results/{model_name}/{sim_name}/{anal_name}/MOSAIC/input/admixedgenofile.22'
	output:
		la_results = "results/{model_name}/{sim_name}/{anal_name}/MOSAIC/localanc_admixed_{nsource}way_1-{naming_mess}.RData",
		model_results = "results/{model_name}/{sim_name}/{anal_name}/MOSAIC/admixed_{nsource}way_1-{naming_mess}.RData",
	log:
		'results/{model_name}/{sim_name}/{anal_name}/MOSAIC/mosaic.{nsource}way_1-{naming_mess}.log'
	params:
		mosaic = config['PATHS']['MOSAIC'],
		input_folder = 'results/{model_name}/{sim_name}/{anal_name}/MOSAIC/input/',
		base_folder = 'results/{model_name}/{sim_name}/{anal_name}/MOSAIC',
		seed = lambda w: units.loc[(w.sim_name, w.anal_name)].anal_seed,
		nthreads = 8,
	shell:
		"""
		Rscript {params.mosaic} --seed {params.seed} --maxcores {params.nthreads} --chromosomes 22:22 --ancestries 3 admixed {params.input_folder} 2>&1 | tee {log}
		# move the results into the results_folder
		# mkdir -p {params.base_folder}/MOSAIC_RESULTS
		# mkdir -p {params.base_folder}/MOSAIC_PLOTS
		# mkdir -p {params.base_folder}/FREQS
		mv MOSAIC_RESULTS/* {params.base_folder}
		mv MOSAIC_PLOTS/* {params.base_folder}
		mv FREQS/* {params.base_folder}
		rmdir MOSAIC_RESULTS
		rmdir MOSAIC_PLOTS
		rmdir FREQS
		"""


rule make_mosaic_input:
	input:
		target_phased_vcf ='results/{model_name}/{sim_name}/{anal_name}/phased.target_inds.vcf.gz',
		reference_vcf ='results/{model_name}/{sim_name}/{anal_name}/phased.reference_inds.vcf.gz',
		plink_map = 'results/{model_name}/{sim_name}/{anal_name}/plink_map.txt',
		site_ts = 'results/{model_name}/{sim_name}/{anal_name}/sample.filter.tsz',
	output:
		'results/{model_name}/{sim_name}/{anal_name}/MOSAIC/input/admixedgenofile.{CHR}'
	params:
		folder = 'results/{model_name}/{sim_name}/{anal_name}/MOSAIC/input',
		chrom_id = lambda w: units.loc[(w.sim_name, w.anal_name)].chr,
		nind_ref = lambda w: units.loc[(w.sim_name, w.anal_name)].nind_ref,
	script:
		'../scripts/make_mosaic_inputs.py'


rule run_RFMix:
	input:
		#alleles_file = 'results/local_ancestry/{model_name}/{sim_name}/{anal_name}.alleles_file.txt',
		alleles_file = 'results/{model_name}/{sim_name}/{anal_name}/MOSAIC/input/admixedgenofile.22',
		classes_file = 'results/{model_name}/{sim_name}/{anal_name}/classes_file.txt',
		positions_file = 'results/{model_name}/{sim_name}/{anal_name}/positions_file.txt',
	output:
		viterbi = 'results/{model_name}/{sim_name}/{anal_name}/RFMix/0.Viterbi.txt'
	params:
		rfmix = config['PATHS']['RFMix'],
		prefix = 'results/{model_name}/{sim_name}/{anal_name}/RFMix/RFMix',
		n_threads = 8,
	shell:
		"""
		cd programs/RFmix/RFMix_v1.5.4
		"""
		"python2 RunRFMix.py TrioPhased "
		"--num-threads {params.n_threads} "
		"../../../{input.alleles_file} ../../../{input.classes_file} ../../../{input.positions_file} "
		"-o ../../../{params.prefix} "


rule make_rfmix_input:
	input:
		site_ts = 'results/{model_name}/{sim_name}/{anal_name}/sample.filter.tsz',
		genetic_map = 'results/{model_name}/{sim_name}/{anal_name}/genetic_map.txt',
	output:
		alleles_file = 'results/{model_name}/{sim_name}/{anal_name}/alleles_file.txt',
		classes_file = 'results/{model_name}/{sim_name}/{anal_name}/classes_file.txt',
		positions_file = 'results/{model_name}/{sim_name}/{anal_name}/positions_file.txt',
	params:
		nind_ref = lambda w: units.loc[(w.sim_name, w.anal_name)].nind_ref,
		nind_admixed = lambda w: units.loc[(w.sim_name, w.anal_name)].nind_admixed,
	script:
		'../scripts/make_rfmix_input.py'


rule run_RFMix2:
	input:
		target_vcf = 'results/{model_name}/{sim_name}/{anal_name}/phased.target_inds.vcf.gz',
		reference_vcf = 'results/{model_name}/{sim_name}/{anal_name}/phased.reference_inds.vcf.gz',
		sample_map = 'results/{model_name}/{sim_name}/{anal_name}/sample_map.txt',
		genetic_map = 'results/{model_name}/{sim_name}/{anal_name}/genetic_map.txt',
	output:
		fb = 'results/{model_name}/{sim_name}/{anal_name}/RFMix2/rfmix2.fb.tsv',
	log:
		'results/{model_name}/{sim_name}/{anal_name}/RFMix2/rfmix2.log',
	params:
		rfmix2 = config['PATHS']['RFMix2'],
		bcftools = config['PATHS']['BCFTOOLS'],
		output = 'results/{model_name}/{sim_name}/{anal_name}/RFMix2/rfmix2',
		nthreads = 8,
		chr = lambda w: simulations.loc[w.sim_name].chr,
		seed = lambda w: units.loc[(w.sim_name, w.anal_name)].anal_seed,
	shell:
		"""
		export PATH="{params.bcftools}:$PATH"
		"""
		"{params.rfmix2} -f {input.target_vcf} -r {input.reference_vcf} "
		"-m {input.sample_map} -g {input.genetic_map} -o {params.output} "
		"--reanalyze-reference -e 5 "
		"--n-threads={params.nthreads} --chromosome={params.chr} "
		"--random-seed={params.seed} 2>&1 | tee {log}"
