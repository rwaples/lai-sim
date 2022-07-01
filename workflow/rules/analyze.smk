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
		# now hardcoded to be 22
		# chrom_id = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].chr,
		nind_ref = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].nind_ref,
		target_pop = lambda w: units.loc[(w.sim_name, w.asc_name)].target_pop[0],
	log:
		'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/MOSAIC/input/make_input.22.log'
	script:
		'../scripts/make_mosaic_inputs.py'


rule run_flare:
	input:
		target_vcf = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/phased.target_inds.vcf.gz',
		reference_vcf = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/phased.reference_inds.vcf.gz',
		sample_map = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/sample_map.txt',
		genetic_map = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/plink_map.txt',
	output:
		'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/flare/flare.anc.vcf.gz',
	params:
		FLARE = config['PATHS']['FLARE'],
		bcftools = config['PATHS']['BCFTOOLS'],
		prefix = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/flare/flare',
		nthreads = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].nthreads,
		seed = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].anal_seed,
	log:
		'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/flare/run_flare.log',
	benchmark:
		'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/benchmark/run_flare.tsv',
	shell:
		"java -Xmx250g -jar {params.FLARE} "
		"ref={input.reference_vcf} "
		"ref-panel={input.sample_map} "
		"gt={input.target_vcf} "
		"map={input.genetic_map} "
		"out={params.prefix} "
		"probs=true "
		"nthreads={params.nthreads} "
		"seed={params.seed} "
		# "min-maf=0 "

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
		'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/logs/make_rfmix_input.log',
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
		# output files that will be gzipped
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
		"-e 5 "
		"--n-threads={params.nthreads} --chromosome={params.chr} "
		"--random-seed={params.seed} 2>&1 | tee {log} "

		"""
		gzip {params.fb}
		gzip {params.Q}
		gzip {params.msp}
		"""


rule run_beagle:
	input:
		gt = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/genotypes.witherror.vcf.gz',
		map = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/plink_map.txt',
	output:
		vcf = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/phased.vcf.gz',
		index = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/phased.vcf.gz.csi',
	log:
		'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/phased.vcf.beagle.log',
	params:
		prefix = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/phased',
		bcftools = config['PATHS']['BCFTOOLS'],
		BEAGLE = config['PATHS']['BEAGLE'],
		seed = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].anal_seed,
	shell:
		"java -jar {params.BEAGLE} "
		"gt={input.gt} "
		"map={input.map} "
		"seed={params.seed} "
		"out={params.prefix} 2>&1 | tee {log} "

		"""

		{params.bcftools} index {output.vcf}
		"""


rule get_dosage_true:
	input:
		true_la = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/true_local_ancestry.site_matrix.npz',
	output:
		'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/ancestry_dosage.true.npz',
	params:
		nsource = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].nsource,
	script:
		'../scripts/get_dosage.true.py'


rule convert_flare_vcf:
	input:
		flare_la = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/flare/flare.anc.vcf.gz',
	output:
		flare_csv = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/flare/flare.anc.csv',
		flare_sites = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/flare/flare.anc.flare_sites',
	params:
		BCFTOOLS = config['PATHS']['BCFTOOLS'],
		nsource = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].nsource,
	shell:
		"""
		{params.BCFTOOLS} query -f '%CHROM, %POS, [%ANP1, %ANP2,]\\n' {input.flare_la} > {output.flare_csv}

		{params.BCFTOOLS} query -f '%POS\n' {input.flare_la}>  {output.flare_sites}
		"""

rule get_dosage_flare:
	input:
		flare_csv = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/flare/flare.anc.csv',
		flare_sites = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/flare/flare.anc.flare_sites',
	output:
		'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/ancestry_dosage.flare.npz',
	params:
		nsource = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].nsource,
		sites_file = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/site.positions',
		BCFTOOLS = config['PATHS']['BCFTOOLS'],
	script:
		'../scripts/get_dosage.flare.py'


rule get_dosage_mosaic:
	input:
		mosaic_la = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/MOSAIC/la_probs.npz',
	output:
		'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/ancestry_dosage.mosaic.npz',
	params:
		nsource = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].nsource,
	script:
		'../scripts/get_dosage.mosaic.py'


rule get_dosage_rfmix2:
	input:
		rfmix2_la = 'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/RFMix2/rfmix2.fb.tsv.gz',
	output:
		'results/{model_name}/{sim_name}/{asc_name}/{anal_name}/ancestry_dosage.rfmix2.npz',
	params:
		nsource = lambda w: units.loc[(w.sim_name, w.asc_name, w.anal_name)].nsource,
	script:
		'../scripts/get_dosage.rfmix2.py'
