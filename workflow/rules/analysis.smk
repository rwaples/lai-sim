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


rule export_mosaic:
	input:
		la_results = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/MOSAIC/localanc_admixed.RData',
		model_results = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/MOSAIC/admixed.RData',
	output:
		path = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/MOSAIC/la_probs.RData',
		np_path = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/MOSAIC/la_probs.npz',
	params:
		input_dir = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/MOSAIC/input/',
	log:
		input_dir = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/MOSAIC/export_mosaic.log',
	script:
		"../scripts/export_mosaic_results.R"


rule run_mosaic:
	input:
		'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/MOSAIC/input/admixedgenofile.22',
	output:
		la_results = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/MOSAIC/localanc_admixed.RData',
		model_results = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/MOSAIC/admixed.RData',
		emlog1 = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/MOSAIC/admixed_1way.EMlog.out',
		emlog3 = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/MOSAIC/admixed_ALLway.EMlog.out',
	params:
		mosaic = config['PATHS']['MOSAIC'],
		# input_folder = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/MOSAIC/input/',
		input_folder = lambda w, input: os.path.dirname(input[0]) + '/',
		base_folder = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/MOSAIC',
		sites_file = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/site.positions',
		seed = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].anal_seed,
		nsource = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].nsource,
		nthreads = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].nthreads,
	log:
		'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/MOSAIC/run_mosaic.log',
	benchmark:
		'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/benchmark/run_mosaic.tsv',
	shadow:
		'full'
	shell:
		"""
		# recalculate gpcm based on the actual number of sites in the analysis
		nsite=$(sed -n '$=' {params.sites_file})
		numer=0.0012
		gpcm=$(python -c "a=$nsite*$numer;print(int(a))")

		Rscript {params.mosaic} --seed {params.seed} --maxcores {params.nthreads} --GpcM $gpcm --chromosomes 22:22 --ancestries {params.nsource} admixed {params.input_folder} 2>&1 | tee {log}

		"""

		"""

		# MOSAIC provides no ability to specify output folders or file names
		# file names depend on the input data and parameters
		# below I try to enforce a consistent naming scheme.

		# move the results into the MOSAIC folder
		mv MOSAIC_RESULTS/* {params.base_folder}
		mv MOSAIC_PLOTS/* {params.base_folder}
		mv FREQS/* {params.base_folder}
		rmdir MOSAIC_RESULTS
		rmdir MOSAIC_PLOTS
		rmdir FREQS

		# try to fix the naming mess here
		cd {params.base_folder}
		find localanc_admixed_*.RData -exec cp {{}} localanc_admixed.RData \;
		find admixed_*.RData -exec cp {{}} admixed.RData \;
		# copy log files
		find admixed_1way*_EMlog.out -exec cp {{}} admixed_1way.EMlog.out \;
		find admixed_{params.nsource}way*_EMlog.out -exec cp {{}} admixed_ALLway.EMlog.out \;

		"""


rule make_mosaic_input:
	input:
		target_phased_vcf = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/phased.target_inds.vcf.gz',
		reference_vcf = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/phased.reference_inds.vcf.gz',
		plink_map = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/plink_map.txt',
		site_ts = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/sample.filter.tsz',
	output:
		'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/MOSAIC/input/admixedgenofile.22',
		'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/MOSAIC/input/pop_0genofile.22',
		'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/MOSAIC/input/pop_1genofile.22',
	params:
		folder = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/MOSAIC/input',
		# chrom_id = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].chr,
		nind_ref = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].nind_ref,
		target_pop = lambda w: units.loc[(w.sim_name, w.asc_name)].target_pop[0],
	log:
		'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/MOSAIC/input/make_input.22.log'
	script:
		'../scripts/make_mosaic_inputs.py'


rule run_bmix:
	input:
		target_vcf = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/phased.target_inds.vcf.gz',
		reference_vcf = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/phased.reference_inds.vcf.gz',
		sample_map = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/sample_map.txt',
		genetic_map = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/plink_map.txt',
	output:
		'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/bmix/bmix.anc.vcf.gz',
	params:
		BMIX = config['PATHS']['BMIX'],
		bcftools = config['PATHS']['BCFTOOLS'],
		prefix = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/bmix/bmix',
		nthreads = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].nthreads,
		seed = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].anal_seed,
	log:
		'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/bmix/run_bmix.log',
	benchmark:
		'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/benchmark/run_bmix.tsv',
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

		{params.bcftools} index {params.prefix}.anc.vcf.gz
		"""

# notice this is RFMix v1, not currently used
rule run_RFMix:
	input:
		alleles_file = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/MOSAIC/input/admixedgenofile.22',
		classes_file = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/classes_file.txt',
		positions_file = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/positions_file.txt',
	output:
		viterbi = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/RFMix/0.Viterbi.txt',
	params:
		rfmix = config['PATHS']['RFMix'],
		prefix = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/RFMix/RFMix',
		nthreads = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].nthreads,
	benchmark:
		'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/benchmark/run_RFMix.tsv',
	log:
		'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/RFMix/run_RFMix.log',
	shell:
		"""
		cd programs/RFmix/RFMix_v1.5.4
		"""
		"python2 RunRFMix.py TrioPhased "
		"--num-threads {params.nthreads} "
		"../../../{input.alleles_file} ../../../{input.classes_file} ../../../{input.positions_file} "
		"-o ../../../{params.prefix} "


rule make_rfmix_input:
	input:
		site_ts = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/sample.filter.tsz',
		genetic_map = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/genetic_map.txt',
	output:
		alleles_file = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/alleles_file.txt',
		classes_file = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/classes_file.txt',
		positions_file = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/positions_file.txt',
	params:
		nind_ref = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].nind_ref,
		nind_admixed = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].nind_admixed,
	log:
		positions_file = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/logs/make_rfmix_input.log',
	script:
		'../scripts/make_rfmix_input.py'


rule run_RFMix2:
	input:
		target_vcf = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/phased.target_inds.vcf.gz',
		reference_vcf = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/phased.reference_inds.vcf.gz',
		sample_map = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/sample_map.txt',
		genetic_map = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/genetic_map.txt',
	output:
		fb = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/RFMix2/rfmix2.fb.tsv.gz',
		Q = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/RFMix2/rfmix2.rfmix.Q.gz',
		msp = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/RFMix2/rfmix2.msp.tsv.gz',
		sis = temp('results/{model_name}/{sim_name}/{asc_name}/{anal_name}/RFMix2/rfmix2.sis.tsv'),
	params:
		rfmix2 = config['PATHS']['RFMix2'],
		bcftools = config['PATHS']['BCFTOOLS'],
		output = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/RFMix2/rfmix2',
		nthreads = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].nthreads,
		chr = lambda w: simulations.loc[w.sim_name].chr,
		seed = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].anal_seed,
		# output files that will be zipped
		fb = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/RFMix2/rfmix2.fb.tsv',
		Q = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/RFMix2/rfmix2.rfmix.Q',
		msp = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/RFMix2/rfmix2.msp.tsv',
	log:
		'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/logs/run_RFMix2.log',
	benchmark:
		'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/benchmark/run_RFMix2.tsv',
	shell:
		"""
		export PATH="{params.bcftools}:$PATH"
		"""
		"{params.rfmix2} -f {input.target_vcf} -r {input.reference_vcf} "
		"-m {input.sample_map} -g {input.genetic_map} -o {params.output} "
		# "--reanalyze-reference "
		# "-e 5 "
		"--n-threads={params.nthreads} --chromosome={params.chr} "
		"--random-seed={params.seed} 2>&1 | tee {log} "

		"""
		gzip {params.fb}
		gzip {params.Q}
		gzip {params.msp}
		"""
