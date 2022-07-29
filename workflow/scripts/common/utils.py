"""General utilty functions used in the Snakemake pipeline."""
import numpy as np
import pandas as pd
import collections
import itertools
import allel
import os
import multiprocessing
import matplotlib.pyplot as plt
import seaborn as sns
import numba


def sample_inds(ts, pop_id, nind, seed):
	"""Return haploid ids based on selecting `nind` diploids from pop_id in ts."""
	rng = np.random.default_rng(seed)

	hap_samples = np.intersect1d(
		ts.samples(population=pop_id),
		np.where(ts.tables.nodes.asdict()['time'] == 0)[0]
	)
	# sample from the first haploids of each ind
	take = rng.choice(hap_samples[::2], nind, replace=False)
	samples = np.empty(nind * 2, dtype=int)
	samples[0::2] = take
	samples[1::2] = take + 1
	return(samples)


def strip_MAC(ts, MAC):
	"""Remove sites with minor allele count <= MAC.

	Returns a new tree-sequence with sites removed
	"""
	initial_sites = ts.num_sites
	present_samples = np.intersect1d(
		ts.samples(),
		np.where(ts.tables.nodes.asdict()['time'] == 0)[0]
	)
	npresent_samples = len(present_samples)

	sites_to_remove = []
	for tree in ts.trees(tracked_samples=present_samples):
		for site in tree.sites():
			assert len(site.mutations) == 1  # Only supports infinite sites muts.
			mut = site.mutations[0]
			if (tree.num_tracked_samples(mut.node) <= MAC) or \
				(tree.num_tracked_samples(mut.node) >= (npresent_samples - MAC)):
				sites_to_remove.append(site.id)

	ts = ts.delete_sites(sites_to_remove)
	final_sites = ts.num_sites
	print(f"MAC filter (<={MAC}):")
	print(f"removed {initial_sites-final_sites} sites ({(initial_sites-final_sites) / (initial_sites):.0%}), {final_sites} sites remain")
	return(ts)


def strip_adjacent_sites(ts, dist=1.5):
	"""Remove sites within dist bp of each other.

	Removes the right-most site in each pair.
	Returns a new tree-sequence.
	"""
	initial_sites = ts.num_sites
	present_samples = np.intersect1d(
		ts.samples(),
		np.where(ts.tables.nodes.asdict()['time'] == 0)[0]
	)
	sites_to_remove = []
	prev_pos = 0
	for tree in ts.trees(tracked_samples=present_samples):
		for site in tree.sites():
			assert len(site.mutations) == 1  # Only supports infinite sites muts.
			# mut = site.mutations[0]
			pos = site.position
			if (pos - prev_pos) < dist:
				sites_to_remove.append(site.id)
			prev_pos = pos

	ts = ts.delete_sites(sites_to_remove)
	final_sites = ts.num_sites
	print('Adjacent sites filter:')
	print(f"removed {initial_sites-final_sites} sites ({(initial_sites-final_sites) / (initial_sites):.0%}), {final_sites} sites remain")
	return(ts)


def downsample_snps(ts, nsnps, seed, fail=False):
	"""Downsample a ts to a maximum number of snps.

	If fail is True, will fail if there are not at least nsnps
	"""
	initial_sites = ts.num_sites
	if ts.num_mutations == nsnps:
		return(ts)
	elif ts.num_mutations < nsnps:
		if fail:
			assert False, "less than {nsnps} are present and fail=True"
		else:
			return(ts)
	else:
		rng = np.random.default_rng(seed)
		keep = frozenset(rng.choice(a=ts.num_sites, size=nsnps, replace=False))
		sites_to_remove = list(frozenset(np.arange(ts.num_sites, dtype=np.int32)) - keep)
		ts = ts.delete_sites(sites_to_remove)
		final_sites = ts.num_sites
		print('Downsample filter:')
		print(f"removed {initial_sites-final_sites} sites ({(initial_sites-final_sites) / (initial_sites):.0%}), {final_sites} sites remain")
		return(ts)


def get_local_ancestry(ts, admixture_time, per_batch):
	"""Return a df describing local ancestry.

	local ancestry is defined by the location of ancestors at admixture_time
	reports local ancestry for nsample chromosomes per population
	and does this by considering per_rep samples at once
	"""
	# target_samples are at time==0
	target_samples = np.intersect1d(
		ts.samples(),
		np.where(ts.tables.nodes.asdict()['time'] == 0)[0]
	)
	# target ancestors exist at time==admixture_time
	tarstors = np.where(ts.tables.nodes.asdict()['time'] == admixture_time)[0]

	nsample = len(target_samples)
	L = [x for x in range(0, nsample, per_batch)]
	R = [x for x in range(per_batch, nsample + per_batch, per_batch)]
	dfs = []
	for i in range(len(L)):
		local = ts.tables.link_ancestors(
			samples=target_samples[L[i]:R[i]],
			ancestors=tarstors
		)

		local_df = pd.DataFrame({
			'left': local.left,
			'right': local.right,
			'parent': local.parent,
			'child': local.child
		})

		dfs.append(local_df)

	local_ancestry_df = pd.concat(dfs)
	pop_of_node = dict(zip(range(len(ts.tables.nodes)), ts.tables.nodes.population))
	# local ancestry population
	local_ancestry_df['localpop'] = [pop_of_node[x] for x in local_ancestry_df['parent']]
	# sampling population
	local_ancestry_df['samplepop'] = [pop_of_node[x] for x in local_ancestry_df['child']]
	local_ancestry_df = local_ancestry_df.sort_values(['samplepop', 'child', 'left']).reset_index(drop=True)
	return(local_ancestry_df)


def get_local_ancestry_pop(ts, pop, admixture_time, max_samples=100, per_rep=12):
	"""Return df describing local ancestry.

	Local ancestry is defined by the location of ancestors at admixture_time
	and does this by considering per_rep samples at once.
	"""
	ancestors = np.where(ts.tables.nodes.asdict()['time'] == admixture_time)[0]
	pop_samples = ts.samples(population=pop)
	if max_samples:
		nsample = np.min([max_samples, pop_samples.size])  # number of samples to report
	else:
		nsample = pop_samples.size
	target_samples = pop_samples[:nsample]  # could allow random sampling of inds

	L = [x for x in range(0, nsample, per_rep)]
	R = [x for x in range(per_rep, nsample, per_rep)] + [nsample]
	assert(len(L) == len(R))
	dfs = []

	for i in range(len(L)):
		local = ts.tables.link_ancestors(
			samples=target_samples[L[i]:R[i]],
			ancestors=ancestors
		)
		local_df = pd.DataFrame({
			'left': local.left,
			'right': local.right,
			'parent': local.parent,
			'child': local.child
		})
		dfs.append(local_df)
		if i % 100 == 0:
			print(f'Done with local ancestry batch {i} out of {len(L)}!')

	local_ancestry_df = pd.concat(dfs)
	pop_of_node = dict(zip(range(len(ts.tables.nodes)), ts.tables.nodes.population))
	# record the population of the ancestor, ie the local ancestry (localpop)
	local_ancestry_df['localpop'] = [pop_of_node[x] for x in local_ancestry_df['parent']]
	# record the population of the sample (samplepop)
	local_ancestry_df['samplepop'] = [pop_of_node[x] for x in local_ancestry_df['child']]
	local_ancestry_df = local_ancestry_df.sort_values(['samplepop', 'child', 'left']).reset_index(drop=True)
	return(local_ancestry_df)


def get_la_mat(ts, df, mapping=None):
	"""Get a array of local ancestry at each site."""
	# the df gives the local ancestry
	# the ts should be simplified with with only the desired (target) samples retained.
	# The node re-mapping produced by simplify(map_nodes=True) should be supplied here if
	# the df was generated prior to simplification.

	# make sure the local ancestry df is sorted by samplepop
	df = df.sort_values(['samplepop', 'child', 'left']).reset_index(drop=True)
	sample_order = df['child'].unique()  # preserves order of occurance
	# this will be the output order of the columns in the la_mat

	# create an empty output matrix and fill with negative ones
	# same size as the genotype matrix
	la_mat = np.zeros((ts.num_sites, ts.num_samples)).astype(int) - 1

	# bed_mask records which sites fall within
	# the intervals on each row of df
	site_pos = np.array([site.position for site in ts.sites()])[:, None]
	bed_mask = (df['left'].values <= site_pos) & (df['right'].values > site_pos)

	sites_idx = np.where(bed_mask)[0]  # gives the site index
	rows_idx = np.where(bed_mask)[1]  # gives the LA row index
	old_of_index = df['child'].to_dict()

	if mapping is not None:
		# use the mapping to get "old" node ids
		new_node, old_node = np.where(mapping == ts.samples()[:, None])
		new_of_old = dict(zip(old_node, new_node))
		# old_of_new = dict(zip(new_node, old_node))
		haps_idx = np.array([new_of_old[old_of_index[x]] for x in rows_idx])  # hap

	else:
		haps_idx = np.array([old_of_index[x] for x in rows_idx])  # hap

	pops = np.array(df['localpop'])[rows_idx]  # pop at each location

	la_mat[sites_idx, haps_idx] = pops
	# finally sort the output by samplepop so the samples from each pop are contiguous
	if mapping is not None:
		poporder = np.array([new_of_old[x] for x in sample_order])
	else:
		poporder = sample_order
	la_mat = la_mat[:, poporder]

	return(la_mat, sample_order)


def get_la_mat_large(ts, df, mapping):
	"""Return a matrix with the local ancestry of each individual at each site.

	the df gives the local ancestry
	the ts should be simplified with with only the target samples retained.
	The node re-mapping produced by simplify(map_nodes=True) should be supplied here.
	"""
	# make sure the local ancestry df is sorted by samplepop
	df = df.sort_values(['samplepop', 'child', 'left']).reset_index(drop=True)
	sample_order = df['child'].unique()  # preserves order of occurance
	# this will be the output order of the columns in the la_mat

	old_of_index = df['child'].to_dict()

	iter_nodes = zip(
		range(len(mapping)),
		mapping
	)
	new_of_old_iter = itertools.filterfalse(lambda x: x[1] == -1, iter_nodes)
	new_of_old = dict(new_of_old_iter)

	pops = df['localpop'].values

	# create an empty output matrix
	# same size as the genotype matrix
	# fill with negative ones
	la_mat = np.zeros((ts.num_sites, ts.num_samples), dtype=np.int8) - 1
	site_pos = np.array([site.position for site in ts.sites()])
	lefts = df['left'].values
	rights = df['right'].values

	for i, pos in enumerate(site_pos):
		sites_idx = i
		site_mask = (lefts <= pos) & (rights > pos)
		rows_idx = np.where(site_mask)[0]
		try:
			haps_idx = np.array([new_of_old[old_of_index[x]] for x in rows_idx])
		except KeyError:
			print(i, pos, rows_idx)
			return(rows_idx)
			haps_idx = np.array([new_of_old[old_of_index[x]] for x in rows_idx])
		site_pops = pops[rows_idx]
		# fill out the la_mat
		la_mat[sites_idx, haps_idx] = site_pops
		if i % 10000 == 0:
			print(f"Done with {i} sites!")

	# finally sort the output by samplepop so the samples from each pop are contiguous
	poporder = np.array([new_of_old[x] for x in sample_order])
	la_mat = la_mat[:, poporder]
	return(la_mat, sample_order)


def find_la(params):
	"""Run link_ancestors() on each set of samples.

	This function to be called from get_local_ancestry_pop_multi()
	"""
	# unpack params
	tables, samples, ancestors = params

	local = tables.link_ancestors(
		samples=samples,
		ancestors=ancestors
	)
	local_df = pd.DataFrame({
		'left': local.left,
		'right': local.right,
		'parent': local.parent,
		'child': local.child
	})
	return(local_df)


def get_local_ancestry_pop_multi(ts, pop, admixture_time, max_samples=100, per_rep=12, n_cores=8):
	"""Return a df describing local ancestry.

	Local ancestry is defined by the location of ancestors at admixture_time
	reports local ancestry for nsample chromosomes per population
	and does this by considering per_rep samples at once.
	"""
	ancestors = np.where(ts.tables.nodes.asdict()['time'] == admixture_time)[0]
	pop_samples = ts.samples(population=pop)
	if max_samples:
		nsample = np.min([max_samples, len(pop_samples)])  # number of samples to report
	else:
		nsample = len(pop_samples)
	print(f'getting local ancestry tracts for {nsample} samples from population {pop}.')
	target_samples = pop_samples[:nsample]  # could allow random sampling of inds

	L = [x for x in range(0, nsample, per_rep)]
	R = [x for x in range(per_rep, nsample, per_rep)] + [nsample]
	assert(len(L) == len(R))

	# multiprocessing starts here
	dfs = []
	pool = multiprocessing.Pool(n_cores)
	table_iter = itertools.repeat(ts.tables.copy(), len(L))
	ancestors_iter = itertools.repeat(ancestors, len(L))

	inputs = zip(
		table_iter,
		[target_samples[L[i]:R[i]] for i in range(len(L))],
		ancestors_iter
	)

	for res in pool.imap(find_la, inputs):
		dfs.append(res)
	pool.close()
	pool.join()

	local_ancestry_df = pd.concat(dfs)
	pop_of_node = dict(zip(range(len(ts.tables.nodes)), ts.tables.nodes.population))
	# record the population of the ancestor, ie the local ancestry (localpop)
	# record the population of the sample (samplepop)
	local_ancestry_df['localpop'] = [pop_of_node[x] for x in local_ancestry_df['parent']]
	# sampling population
	local_ancestry_df['samplepop'] = [pop_of_node[x] for x in local_ancestry_df['child']]
	local_ancestry_df = local_ancestry_df.sort_values(['samplepop', 'child', 'left']).reset_index(drop=True)
	return(local_ancestry_df)


def make_ind_labels(ts):
	"""Make labels for diploid samples present in the ts."""
	pop_of_sample = dict(zip(range(len(ts.tables.nodes)), ts.tables.nodes.population))
	pops = [pop_of_sample[i] for i in ts.samples()]
	ind_labels = []
	countpop = collections.defaultdict(int)
	for indpop in pops[::2]:
		countpop[indpop] += 1
		label = f'pop_{indpop}-ind_{countpop[indpop]:04}'
		ind_labels.append(label)
	return(ind_labels)


def vcfheader(ts, target_pop):
	"""Construct an appropriate vcf header.

	Always uses chr22 as contig ID.
	Includes current date.
	"""
	from datetime import date
	fileformat = '##fileformat=VCFv4.2\n'
	fileDate = f'##filedate={date.today().ctime()}\n'
	FILTER = '##FILTER=<ID=PASS,Description="All filters passed">\n'
	contig = f'##contig=<ID=chr22,length={int(np.ceil(ts.get_sequence_length()))}>\n'
	FORMAT = '##FORMAT=<ID=LA,Number=1,Type=String,Description="Local ancestry">\n'

	ind_labels = make_ind_labels(ts.simplify(ts.samples(target_pop)))

	header = fileformat + fileDate + FILTER + contig + FORMAT
	header = header + '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'])
	header = header + '\t' + '\t'.join(ind_labels) + '\n'
	return(header)


def get_allele_freqs(ts, pops=None):
	"""Return allele frequencies for each site and population.

	returns a matrix of [site_idx, pop_idx] derived allele frequencies.
	pops is a list of the pops, order of pops is maintained.
	Only supports infinite site mutations with <=2 alleles.
	"""
	if pops:
		freqs = np.zeros((ts.num_sites, len(pops)))
	else:
		pops = [p.id for p in ts.populations()]
		freqs = np.zeros((ts.num_sites, len(pops)))
	for ipop, pop in enumerate(pops):
		samples = ts.samples(pop)
		nsamp = len(samples)
		isite = 0
		for tree in ts.trees(tracked_samples=samples):
			for site in tree.sites():
				assert len(site.mutations) == 1  # Only supports infinite sites
				mut = site.mutations[0]
				ac = tree.num_tracked_samples(mut.node)
				af = ac / nsamp
				freqs[isite, ipop] = af
				isite += 1
	return(freqs)


def get_mean_allele_frequencies(freqs):
	"""Return the mean allele frquency across populations.

	Input should be the array from get_allele_freqs().
	"""
	return(freqs.mean(1))


def get_fst_faser(ts, popA, popB):
	"""Return Fst (Hudson) for popA and popB.

	Uses the allel.hudson_fst() function from scikit-allel.
	popA and popB give the samples to use from each pop.
	"""
	popA_samples = ts.samples(popA)
	popB_samples = ts.samples(popB)
	# count the number of derived alleles within each pop
	# without the 'int64' conversion here
	# there seems to be an overflow error when calculating Fst with many samples
	popA_counts = np.zeros(ts.num_sites, dtype='int64')
	popB_counts = np.zeros(ts.num_sites, dtype='int64')
	s = 0
	for tree in ts.trees(tracked_samples=popA_samples):
		for site in tree.sites():
			assert len(site.mutations) == 1  # Only supports infinite sites muts.
			mut = site.mutations[0]
			popA_counts[s] = tree.num_tracked_samples(mut.node)
			s += 1
	s = 0
	for tree in ts.trees(tracked_samples=popB_samples):
		for site in tree.sites():
			assert len(site.mutations) == 1  # Only supports infinite sites muts.
			mut = site.mutations[0]
			popB_counts[s] = tree.num_tracked_samples(mut.node)
			s += 1

	npopA = len(popA_samples)
	npopB = len(popB_samples)
	# get the number of ancestral and derived alleles within each pop
	acA = np.array([npopA - popA_counts, popA_counts]).T
	acB = np.array([npopB - popB_counts, popB_counts]).T
	num, denom = allel.hudson_fst(acA, acB)
	fst = np.sum(num) / np.sum(denom)
	return(fst)


def get_ancestry_dosage(arr, n_anc):
	"""Compute ancestry dosage from probablistic haploid ancestry calls."""
	#anc_dosage = np.zeros((arr.shape[0], int(arr.shape[1] / 2)), dtype=np.half)
	anc_dosage = np.zeros((arr.shape[0], int(arr.shape[1] / 2)))
	if n_anc == 2:
		a0 = arr[:, 0::2]  # should be views
		a1 = arr[:, 1::2]
		anc_dosage[:, 0::2] = a0[:, ::2] + a0[:, 1::2]
		anc_dosage[:, 1::2] = a1[:, ::2] + a1[:, 1::2]
	if n_anc == 3:
		a0 = arr[:, 0::3]
		a1 = arr[:, 1::3]
		a2 = arr[:, 2::3]
		anc_dosage[:, 0::3] = a0[:, ::2] + a0[:, 1::2]
		anc_dosage[:, 1::3] = a1[:, ::2] + a1[:, 1::2]
		anc_dosage[:, 2::3] = a2[:, ::2] + a2[:, 1::2]
	elif n_anc == 4:
		a0 = arr[:, 0::4]
		a1 = arr[:, 1::4]
		a2 = arr[:, 2::4]
		a3 = arr[:, 3::4]
		anc_dosage[:, 0::4] = a0[:, ::2] + a0[:, 1::2]
		anc_dosage[:, 1::4] = a1[:, ::2] + a1[:, 1::2]
		anc_dosage[:, 2::4] = a2[:, ::2] + a2[:, 1::2]
		anc_dosage[:, 3::4] = a3[:, ::2] + a3[:, 1::2]
	return(anc_dosage)


def load_true_la(path):
	"""Load true local ancestry."""
	return(np.load(path)['arr'])
	# return(np.load(path)['arr'].astype(np.single))


def get_true_anc_dosage(true_la, n_anc):
	"""Get true local ancestry dosage."""
	hap1 = np.zeros((true_la.shape[0], int(true_la.shape[1] / 2 * n_anc)), dtype='int8')
	hap2 = np.zeros((true_la.shape[0], int(true_la.shape[1] / 2 * n_anc)), dtype='int8')
	aa = np.arange(true_la[:, ::2].shape[1]) * n_anc + true_la[:, ::2]
	bb = np.arange(true_la[:, 1::2].shape[1]) * n_anc + true_la[:, 1::2]
	np.put_along_axis(hap1, aa.astype('int'), 1, axis=1)
	np.put_along_axis(hap2, bb.astype('int'), 1, axis=1)
	return(hap1 + hap2)


@numba.jit(
	numba.float64(
		numba.float32[:],
		numba.float32[:],
	),
	nopython=True
)
def pearsonr2_numba(x, y):
	"""Return the *squared* pearson correlation coef.

	x and y are not modified.

	Input arrays should be np.float32, output is a single float64 value.
	"""
	assert len(x) == len(y)
	# return Nan if either x or y is constant
	if (x == x[0]).all() or (y == y[0]).all():
		return(np.nan)

	n = len(x)
	# mean of x and y
	xmean = numba.float64(0.0)
	ymean = numba.float64(0.0)
	for i in range(n):
		xmean += x[i]
		ymean += y[i]
	xmean = xmean / n
	ymean = ymean / n

	r_num = numba.float64(0.0)
	r_dena = numba.float64(0.0)
	r_denb = numba.float64(0.0)
	for i in range(n):
		xd = x[i] - xmean  # difference from the mean
		yd = y[i] - ymean
		r_num += xd * yd
		r_dena += xd * xd
		r_denb += yd * yd
	r2 = (r_num * r_num) / (r_dena * r_denb)
	return(r2)


def r2_dosage_ancestry(true_dosage, pred_dosage, n_anc):
	"""Get ancestry-specific R2 values for LA vs truth."""
	per_anc = []
	for i in range(n_anc):
		per_anc.append(
			pearsonr2_numba(
				true_dosage[:, i::n_anc].reshape(-1),
				pred_dosage[:, i::n_anc].reshape(-1),
			)
		)
	return(per_anc)


def r2_dosage_individual(true_dosage, pred_dosage, n_anc):
	"""Get individual-specific R2 values for LA vs truth."""
	per_ind = []
	for i in range(int(true_dosage.shape[1] / n_anc)):
		per_ind.append(
			pearsonr2_numba(
				true_dosage[:, i * n_anc:i * n_anc + n_anc].reshape(-1),
				pred_dosage[:, i * n_anc:i * n_anc + n_anc].reshape(-1),
			)
		)
	return(per_ind)


def load_rfmix_fb(path, sites_file):
	"""Load an array of the posterior local ancestry probabilities from RFMixv2."""
	rfmix_res = pd.read_csv(path, sep='\t', comment='#')
	# expand out to each site
	# needed because RFMix2 only reports LA every fifth site.
	rfmix_res = np.repeat(rfmix_res.iloc[:, 4:].values, [5], axis=0)
	# rfmix_res = rfmix_res.astype(np.half)
	# rfmix_res = rfmix_res.astype(np.single)

	nsites = len(pd.read_csv(sites_file, header=None))
	rfmix_res = rfmix_res[:nsites, :]
	return(rfmix_res)


def load_flare(path, sites_file, flare_sites, BCFTOOLS):
	"""Load an array of the posterior local ancestry probabilities from flare."""
	flare = pd.read_csv(path, header=None)
	flare = flare.dropna(axis=1)
	res = flare.iloc[:, 2:].values
	res = np.concatenate([res[:1], res])
	# res = res.astype(np.half)
	# res = res.astype(np.single)

	# account for any sites filtered by flare (e.g. due to MAF)
	pre_sites = pd.read_csv(sites_file, header=None).values.flatten()
	post_sites = pd.read_csv(flare_sites, header=None).values.flatten()
	post_indexes = np.searchsorted(post_sites, pre_sites)
	del pre_sites, post_sites
	res = res[post_indexes]

	return(res)


def load_mosaic(path):
	"""Return an array of the posterior LA probabilities from MOSAIC."""
	arr = np.load(path)['arr']
	res = arr.T.reshape((arr.shape[2], -1), order='C')
	# res = res.astype(np.half)
	# res = res.astype(np.single)
	return(res)


def get_Q(arr, n_anc):
	"""
	Return a data frame of `global` ancestry fractions (Q).

	Calculated as sums over probabalistic local ancestry proportions.
	"""
	nsites = arr.shape[0]
	# avoid overflow and sum over sites
	arr = arr.astype(float).sum(0)
	if n_anc == 2:
		a0 = arr[0::2]  # should be views
		a1 = arr[1::2]
		q0 = a0 / (nsites * 2)
		q1 = a1 / (nsites * 2)
		Q = pd.DataFrame([q0, q1]).T
		Q.columns = ['pop_0', 'pop_1']
	elif n_anc == 3:
		a0 = arr[0::3]
		a1 = arr[1::3]
		a2 = arr[2::3]
		q0 = a0 / (nsites * 2)
		q1 = a1 / (nsites * 2)
		q2 = a2 / (nsites * 2)
		Q = pd.DataFrame([q0, q1, q2]).T
		Q.columns = ['pop_0', 'pop_1', 'pop_2']
	elif n_anc == 4:
		a0 = arr[0::4]
		a1 = arr[1::4]
		a2 = arr[2::4]
		a3 = arr[3::4]
		q0 = a0 / (nsites * 2)
		q1 = a1 / (nsites * 2)
		q2 = a2 / (nsites * 2)
		q3 = a3 / (nsites * 2)
		Q = pd.DataFrame([q0, q1, q2, q3]).T
		Q.columns = ['pop_0', 'pop_1', 'pop_2', 'pop_3']

	return(Q)


def get_RMSD_Q(Q1, Q2):
	"""Return the RMSD between two sets of ancestry proportions (Q).

	RMSD = root mean squared deviation.
	"""
	assert(Q1.shape == Q2.shape)
	D = Q1 - Q2
	SD = D * D
	MSD = SD.mean().mean()
	RMSD = np.sqrt(MSD)
	return(RMSD)


def max_la(vals, n_anc):
	"""
	Replace each probabalistic local ancestry call with a `hard` call.

	Calls the ancestry with the highest posterior probability.
	All other ancestries are assigned probability of zero.
	Breaks ties by going to the lower-numbered population.
	NOTICE - Modifies the passed array in place.
	"""
	idxs = np.arange(0, vals.shape[1] + n_anc, n_anc)
	for i in range(1, len(idxs)):
		b = vals[:, idxs[i - 1]:idxs[i]]
		dims = list(b.shape)
		dims[1] = 1
		a = b.argmax(1).reshape(dims)  # breaks ties by assigning the lower pop
		c = np.zeros_like(b)
		np.put_along_axis(c, a, 1, axis=1)
		vals[:, idxs[i - 1]:idxs[i]] = c
	return(vals)


def plot_ancestry_dosage(
	pred_dosage,
	start_index,
	n_anc, title,
	path=None,
	format='pdf',
	reference_dosage=None
):
	"""Plot ancestry dosage."""
	colors = ['blue', 'orange', 'green', 'purple']

	fig, ax = plt.subplots(
		figsize=(12, n_anc * 1.5),
		nrows=n_anc,
		sharex=True,
		sharey=True
	)
	f = []
	for i in range(n_anc):
		l, = ax[i].plot(pred_dosage[:, start_index + i], c=colors[i])
		f.append(l)
		if reference_dosage is not None:
			l, = ax[i].plot(
				reference_dosage[:, start_index + i], c=colors[i], alpha=.3, ls='--'
			)

	plt.legend(f, [f'pop{p}' for p in range(n_anc)])

	fig.tight_layout()
	sns.despine(bottom=True)
	ax[0].set_title(title)
	ax[-1].set_xlabel('Site number ')
	if path:
		plt.savefig(path, dpi=300, format=format, bbox_inches='tight')


def make_qq_report(inferred_dosage, true_dosage, nbins):
	"""Generate a report of observed vs inferred mean ancestry dosage.

	Bins are equally spaced between 0 to 2.
	"""
	MAXSIZE = 1800000000

	try:
		assert true_dosage.shape == inferred_dosage.shape
	except AssertionError:
		report = pd.DataFrame({'bin': range(nbins), 'mean': 0, 'n': 0})
		return(report)
	qq = pd.DataFrame(
		{
			'true': true_dosage.flatten()[:MAXSIZE],
			'inferred': inferred_dosage.flatten()[:MAXSIZE]
		},
		dtype="float64"
	)
	del inferred_dosage
	del true_dosage

	qq['bin'] = pd.cut(
		qq['inferred'],
		#bins=np.linspace(0, 2, nbins + 1, dtype="float32"),
		bins=np.linspace(0, 2, nbins + 1, dtype="float64"),
		include_lowest=True,
		labels=False
	)
	report = qq.groupby(['bin'])['true'].agg([np.mean, len]).reset_index()
	report.columns = ['bin', 'mean', 'n']
	return(report)


def add_reports(report_a, report_b):
	"""Add two QQ reports and return a new report.

	New report is a weighted average of the two inputs.
	"""

	merged = report_a.merge(report_b, on='bin', how='outer')
	merged = merged.replace(np.nan, 0)
	# n is a sum of a and b
	merged['n'] = np.nansum([merged['n_x'], merged['n_y']], axis=0, dtype=int)
	# mean is a weighted average of a and b
	merged['mean'] = np.average(
		[merged['mean_x'], merged['mean_y']],
		axis=0,
		weights=[merged['n_x'], merged['n_y']]
	)

	return(merged[['bin', 'mean', 'n']])


def plot_qq_report(report, title=None, reflines=True):
	"""Plot a QQ report."""
	if reflines:
		point1 = [-.1, 0]
		point2 = [.1, 0]
		x_values = [point1[0], point2[0]]
		y_values = [point1[1], point2[1]]
		plt.plot(x_values, y_values, c='k', ls='--')
		point1 = [.9, 1]
		point2 = [1.1, 1]
		x_values = [point1[0], point2[0]]
		y_values = [point1[1], point2[1]]
		plt.plot(x_values, y_values, c='k', ls='--')
		point1 = [1.9, 2]
		point2 = [2.1, 2]
		x_values = [point1[0], point2[0]]
		y_values = [point1[1], point2[1]]
		plt.plot(x_values, y_values, c='k', ls='--')
	plt.scatter(
		report['bin'] / len(report) * 2,
		report['mean'],
		c=['r'] + ['b'] * int(np.floor((len(report) - 3) / 2)) + ['r'] + ['b'] * int(np.ceil((len(report) - 3) / 2)) + ['r']
	)
	plt.title(title)
	plt.xlabel('inferred ancestry dosage bin')
	plt.ylabel('mean of true ancestry dosage')
	plt.show()
	plt.scatter(
		report['bin'] / len(report) * 2,
		report['n'],
		c=['r'] + ['b'] * int(np.floor((len(report) - 3) / 2)) + ['r'] + ['b'] * int(np.ceil((len(report) - 3) / 2)) + ['r']
	)

	plt.xlabel('inferred ancestry dosage bin')
	plt.ylabel('count')
	plt.semilogy()
	plt.show()
