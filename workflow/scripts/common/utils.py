import numpy as np
import pandas as pd
import collections
import itertools
import allel
import pyreadr
import os
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
from sklearn.metrics import r2_score


def sample_inds(ts, pop_id, nind, seed):
	"""return the haploid sample ids representing sampling nind individuals from pop_id in ts"""
	rng = np.random.default_rng(seed)

	hap_samples = np.intersect1d(
		ts.samples(population = pop_id),
		np.where(ts.tables.nodes.asdict()['time']==0)[0]
	)
	# sample from the first haploids of each ind
	take = rng.choice(hap_samples[::2], nind, replace=False)
	samples = np.empty(nind*2, dtype=int)
	samples[0::2] = take
	samples[1::2] = take+1
	return samples


def strip_MAC(ts, MAC):
	"""
	Removes sites with minor allele count <= MAC.
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
	print(f"removed {initial_sites-final_sites} sites ({(initial_sites-final_sites)/(initial_sites):.0%}), {final_sites} sites remain")
	return ts


def strip_adjacent_sites(ts, dist=1.5):
	"""
	Removes sites within dist bp of each other.
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
	print(f"removed {initial_sites-final_sites} sites ({(initial_sites-final_sites)/(initial_sites):.0%}), {final_sites} sites remain")
	return ts


def downsample_snps(ts, nsnps, seed, fail=False):
	"""downsample a tree-sequence to a maximum number of snps
	if fail is True, will fail if there are not at least nsnps"""
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
		keep = frozenset(np.random.choice(a=ts.num_mutations, size=nsnps, replace=False))

		# RKW new way to do this
		keep = frozenset(np.random.choice(a=ts.num_sites, size=nsnps, replace=False))
		sites_to_remove = list(frozenset(np.arange(ts.num_sites, dtype = np.int32)) - keep)
		ts = ts.delete_sites(sites_to_remove)
		final_sites = ts.num_sites
		print('Downsample filter:')
		print(f"removed {initial_sites-final_sites} sites ({(initial_sites-final_sites)/(initial_sites):.0%}), {final_sites} sites remain")
		return(ts)


		tables = ts.dump_tables()
		tables.sites.clear()
		tables.mutations.clear()
		mut_ix = 0
		for tree in ts.trees():
			for site in tree.sites():
				nmut = len(site.mutations)
				if(nmut==0):
					pass
				else:
					assert nmut == 1, f"{site}"  # Only supports infinite sites muts.
					if mut_ix in keep:
						mut = site.mutations[0]
						site_id = tables.sites.add_row(
							position=site.position,
							ancestral_state=site.ancestral_state)
						tables.mutations.add_row(
							site=site_id, node=mut.node, derived_state=mut.derived_state)
					mut_ix +=1
		final_sites = len(tables.sites)
		print('Downsample filter:')
		print(f"removed {initial_sites-final_sites} sites ({(initial_sites-final_sites)/(initial_sites):.0%}), {final_sites} sites remain")
		return tables.tree_sequence()


def get_local_ancestry(ts, admixture_time, per_batch):
	"""returns df describing local ancestry
	local ancestry is defined by the location of ancestors at admixture_time
	reports local ancestry for nsample chromosomes per population
	and does this by considering per_rep samples at once
	"""


	# target_samples are at time==0
	target_samples = np.intersect1d(
				ts.samples(),
				np.where(ts.tables.nodes.asdict()['time']==0)[0]
			)
	# target ancestors are at time==admixture_time
	target_ancestors = np.where(ts.tables.nodes.asdict()['time']==admixture_time)[0]

	nsample = len(target_samples)
	l = [x for x in range(0, nsample, per_batch)]
	r = [x for x in range(per_batch, nsample+per_batch, per_batch)]
	dfs = []
	for i in range(len(l)):
		local = ts.tables.link_ancestors(
			samples = target_samples[l[i]:r[i]],
			ancestors = target_ancestors
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
	return local_ancestry_df


def get_local_ancestry_pop(ts, pop, admixture_time, max_samples=100, per_rep=12):
	"""returns df describing local ancestry
	local ancestry is defined by the location of ancestors at admixture_time
	and does this by considering per_rep samples at once
	"""

	ancestors = np.where(ts.tables.nodes.asdict()['time'] == admixture_time)[0]
	pop_samples = ts.samples(population=pop)
	if max_samples:
		nsample = np.min([max_samples, pop_samples.size])  # number of samples to report
	else:
		nsample = pop_samples.size
	target_samples = pop_samples[:nsample]  # could allow random sampling of inds

	l = [x for x in range(0, nsample, per_rep)]
	r = [x for x in range(per_rep, nsample, per_rep)] + [nsample]
	assert(len(l) == len(r))
	dfs = []

	for i in range(len(l)):
		local = ts.tables.link_ancestors(
			samples=target_samples[l[i]:r[i]],
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
			print(f'Done with local ancestry batch {i} out of {len(l)}!')

	local_ancestry_df = pd.concat(dfs)
	pop_of_node = dict(zip(range(len(ts.tables.nodes)), ts.tables.nodes.population))
	# record the population of the ancestor, ie the local ancestry (localpop)
	# record the population of the sample (samplepop)
	local_ancestry_df['localpop'] = [pop_of_node[x] for x in local_ancestry_df['parent']]
	# sampling population
	local_ancestry_df['samplepop'] = [pop_of_node[x] for x in local_ancestry_df['child']]
	local_ancestry_df = local_ancestry_df.sort_values(['samplepop', 'child', 'left']).reset_index(drop=True)
	return local_ancestry_df


def get_la_mat(ts, df, mapping=None):
	# the df gives the local ancestry
	# the ts should be simplified with with only the desired (target) samples retained.
	# The node re-mapping produced by simplify(map_nodes=True) should be supplied here if
	# the df was generated prior to simplification.

	# make sure the local ancestry df is sorted by samplepop
	df = df.sort_values(['samplepop', 'child', 'left']).reset_index(drop=True)
	sample_order = df['child'].unique() # preserves order of occurance
	# this will be the output order of the columns in the la_mat

	# create an empty output matrix
	# fill with negative ones
	# same size as the genotype matrix
	la_mat = np.zeros((ts.num_sites, ts.num_samples)).astype(int) - 1

	# bed_mask records which sites fall within
	# the intervals on each row of df
	site_pos = np.array([site.position for site in ts.sites()])[:,None]
	bed_mask = (df['left'].values <= site_pos) & (df['right'].values > site_pos)

	sites_idx = np.where(bed_mask)[0] # gives the site index
	rows_idx = np.where(bed_mask)[1] # gives the LA row index
	old_of_index = df['child'].to_dict()

	if mapping is not None:
		# use the mapping to get "old" node ids
		new_node, old_node = np.where(mapping == ts.samples()[:,None])
		new_of_old = dict(zip(old_node, new_node))
		old_of_new = dict(zip(new_node, old_node))
		haps_idx = np.array([new_of_old[old_of_index[x]] for x in rows_idx]) # gives the hap

	else:
		haps_idx = np.array([old_of_index[x] for x in rows_idx]) # gives the hap

	pops = np.array(df['localpop'])[rows_idx] # pop at each location

	la_mat[sites_idx, haps_idx] = pops
	# finally sort the output by samplepop so the samples from each pop are contiguous
	if mapping is not None:
		poporder = np.array([new_of_old[x] for x in sample_order])
	else:
		poporder = sample_order
	la_mat = la_mat[:, poporder]

	return la_mat, sample_order


def get_la_mat_large(ts, df, mapping):
	"""
	Returns a matrix with the local ancestry of each individual at each site.
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

	return la_mat, sample_order

def find_la(params):
	"""
	runs link_ancestors() on each set of samples.
	called from get_local_ancestry_pop_multi()
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
	return local_df


def get_local_ancestry_pop_multi(ts, pop, admixture_time, max_samples=100, per_rep=12, n_cores=8):
	"""returns df describing local ancestry
	local ancestry is defined by the location of ancestors at admixture_time
	reports local ancestry for nsample chromosomes per population
	and does this by considering per_rep samples at once
	"""

	ancestors = np.where(ts.tables.nodes.asdict()['time'] == admixture_time)[0]
	pop_samples = ts.samples(population=pop)
	if max_samples:
		nsample = np.min([max_samples, len(pop_samples)])  # number of samples to report
	else:
		nsample = len(pop_samples)
	print(f'getting local ancestry tracts for {nsample} samples from population {pop}.')
	target_samples = pop_samples[:nsample]  # could allow random sampling of inds

	l = [x for x in range(0, nsample, per_rep)]
	r = [x for x in range(per_rep, nsample, per_rep)] + [nsample]
	assert(len(l) == len(r))

	# multiprocessing starts here
	dfs = []
	pool = multiprocessing.Pool(n_cores)
	table_iter = itertools.repeat(ts.tables.copy(), len(l))
	ancestors_iter = itertools.repeat(ancestors, len(l))

	inputs = zip(
		table_iter,
		[target_samples[l[i]:r[i]] for i in range(len(l))],
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
	pop_of_sample = dict(zip(range(len(ts.tables.nodes)), ts.tables.nodes.population))
	nind_in_pop = collections.defaultdict(int)
	pops = [pop_of_sample[i] for i in ts.samples()]
	for p in pops:
		nind_in_pop[p] +=1
	for p in nind_in_pop:
		nind_in_pop[p] = int(nind_in_pop[p]/2)
	ind_labels = []
	for p in nind_in_pop:
		for i in range(1, nind_in_pop[p]+1):
			ind_labels.append(f'pop_{p}-ind_{i:04}')
	ind_labels = sorted(ind_labels)
	return(ind_labels)


def vcfheader(ts, target_pop):
	from datetime import date
	fileformat = '##fileformat=VCFv4.2\n'
	fileDate = f'##filedate={date.today().ctime()}\n'
	FILTER = '##FILTER=<ID=PASS,Description="All filters passed">\n'
	contig = f'##contig=<ID=chr22,length={int(np.ceil(ts.get_sequence_length()))}>\n'
	FORMAT = '##FORMAT=<ID=LA,Number=1,Type=String,Description="Local ancestry">\n'

	pop_of_sample = dict(zip(range(len(ts.tables.nodes)), ts.tables.nodes.population))

	nind_in_pop = collections.defaultdict(int)
	pops = [pop_of_sample[i] for i in ts.samples(target_pop)]
	for p in pops:
		nind_in_pop[p] += 1
	for p in nind_in_pop:
		nind_in_pop[p] = int(nind_in_pop[p]/2)

	ind_labels = []
	for p in nind_in_pop:
		for i in range(1, nind_in_pop[p]+1):
			ind_labels.append(f'pop_{p}-ind_{i:05}')

	header = fileformat+fileDate+FILTER+contig+FORMAT
	header = header +'\t'.join(['#CHROM','POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'])
	header = header + '\t' + '\t'.join(ind_labels) + '\n'
	return(header)


def get_allele_freqs(ts, pops = None):
	"""return allele frequencies for each site and population
	returns a matrix of [site_idx, pop_idx] derived allele frequencies.
	pops is a list of the pops, order of pops is maintained.
	Only supports infinite site mutations with <=2 alleles.
	"""
	if pops:
		freqs = np.zeros((ts.num_sites, len(pops)))
	else:
		pops = [p.id for p in ts2.populations()]
		freqs = np.zeros((ts.num_sites, len(pops)))
	for ipop, pop in enumerate(pops):
		samples = ts.samples(pop)
		nsamp = len(samples)
		isite = 0
		for tree in ts.trees(tracked_samples=samples): # Only supports infinite sites muts.
			for site in tree.sites():
				assert len(site.mutations) == 1
				mut = site.mutations[0]
				ac = tree.num_tracked_samples(mut.node)
				af = ac/nsamp
				freqs[isite, ipop] = af
				isite+=1
	return(freqs)


def get_mean_allele_frequencies(freqs):
	"""return the mean allele frquency across populations.
	Input should be the array from get_allele_freqs(). """
	return(freqs.mean(1))


def get_fst_faser(ts, popA, popB):
	"""
	Hudson Fst for popA and popB.
	Uses the allel.hudson_fst() function from scikit-allel.
	popA and popB are the samples from each pop
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
	anc_dosage = np.zeros((arr.shape[0], int(arr.shape[1]/2)), dtype=np.half)
	if n_anc==2:
		a0 = arr[:, 0::2] # should be views
		a1 = arr[:, 1::2]
		anc_dosage[:, 0::2] = a0[:, ::2] + a0[:, 1::2]
		anc_dosage[:, 1::2] = a1[:, ::2] + a1[:, 1::2]
	if n_anc==3:
		assert (n_anc==3)
		a0 = arr[:, 0::3] # should be views
		a1 = arr[:, 1::3]
		a2 = arr[:, 2::3]
		anc_dosage[:, 0::3] = a0[:, ::2] + a0[:, 1::2]
		anc_dosage[:, 1::3] = a1[:, ::2] + a1[:, 1::2]
		anc_dosage[:, 2::3] = a2[:, ::2] + a2[:, 1::2]
	elif n_anc==4:
		assert (n_anc==4)
		a0 = arr[:, 0::4] # should be views
		a1 = arr[:, 1::4]
		a2 = arr[:, 2::4]
		a3 = arr[:, 3::4]
		anc_dosage[:, 0::4] = a0[:, ::2] + a0[:, 1::2]
		anc_dosage[:, 1::4] = a1[:, ::2] + a1[:, 1::2]
		anc_dosage[:, 2::4] = a2[:, ::2] + a2[:, 1::2]
		anc_dosage[:, 3::4] = a3[:, ::2] + a3[:, 1::2]
	return(anc_dosage)


def load_true_la(path):
	return(np.load(path)['arr'])


def get_true_anc_dosage(true_la, n_anc):
	hap1 = np.zeros((true_la.shape[0], int(true_la.shape[1]/2*n_anc)), dtype = 'int8')
	hap2 = np.zeros((true_la.shape[0], int(true_la.shape[1]/2*n_anc)), dtype = 'int8')
	aa = np.arange(true_la[:, ::2].shape[1])*n_anc+true_la[:, ::2]
	bb = np.arange(true_la[:, 1::2].shape[1])*n_anc+true_la[:, 1::2]
	np.put_along_axis(hap1, aa, 1, axis=1)
	np.put_along_axis(hap2, bb, 1, axis=1)
	return(hap1+hap2)


def r2_ancestry_dosage(true_dosage, pred_dosage, n_anc):
	per_anc = []
	for i in range(n_anc):
		per_anc.append(
			pearsonr(
				true_dosage[:,i::n_anc].ravel(),
				pred_dosage[:,i::n_anc].ravel()
			)[0]
		)
	per_ind = []
	for i in range(int(true_dosage.shape[1]/n_anc)):
		per_ind.append(
			pearsonr(
				true_dosage[:, i*n_anc:i*n_anc+n_anc].ravel(),
				pred_dosage[:, i*n_anc:i*n_anc+n_anc].ravel()
			)[0]
		)

	return(per_anc, per_ind)


## Load in the probablistic output of each method
def load_rfmix_fb(path):
	"""Load and return an array of the posterior local ancestry probabilities from RFMixv2."""
	rfmix_res = pd.read_csv(path, sep='\t', comment='#')
	# expand out to each site
	rfmix_res = np.repeat(rfmix_res.iloc[:, 4:].values, [5], axis = 0)
	return(rfmix_res)


def load_bmix(path, sites_file, BCFTOOLS):
	"""Load and return an array of the posterior local ancestry probabilities from bmix."""

	# convert the vcf.gz to a csv
	csv_path = path.replace('.vcf.gz', '.csv')
	bmix_sites = path.replace('.vcf.gz', '.bmix_sites')
	os.system(f"{BCFTOOLS} query -f '%CHROM, %POS, [%ANP1, %ANP2,]\\n' {path} > {csv_path}")
	os.system(f"{BCFTOOLS} query -f '%POS\n' {path} > {bmix_sites}")

	bmix = pd.read_csv(csv_path, header=None)
	bmix = bmix.dropna(axis=1)
	res = bmix.iloc[:,2:].values
	res = np.concatenate([res[:1], res])

	pre_sites = pd.read_csv(sites_file, header=None).values.flatten()
	post_sites = pd.read_csv(bmix_sites, header=None).values.flatten()
	post_indexes = np.searchsorted(post_sites, pre_sites)
	res = res[post_indexes]
	return(res)


def load_mosaic(path):
	"""Load and return an array of the posterior local ancestry probabilities from MOSAIC."""
	#mr = pyreadr.read_r(path)['arr'].astype(np.half)
	#res = mr.to_numpy().T.reshape((mr.shape[2],-1), order='C')
	arr = np.load(path)['arr']
	res = arr.T.reshape((arr.shape[2],-1), order='C')
	return(res)


def get_Q(arr, n_anc):
	"""
	Return a data frame of ancestry fractions (Q)
	calculated from probabalistic local ancestry proportions.
	"""
	nsites = arr.shape[0]
	# avoid overflow and sum over sites
	arr = arr.astype(float).sum(0)
	if n_anc == 2:
		a0 = arr[0::2] # should be views
		a1 = arr[1::2]
		q0 = a0/(nsites*2)
		q1 = a1/(nsites*2)
		Q = pd.DataFrame([q0, q1]).T
		Q.columns = ['pop_0', 'pop_1']
	elif n_anc == 3:
		a0 = arr[0::3] # should be views
		a1 = arr[1::3]
		a2 = arr[2::3]
		q0 = a0/(nsites*2)
		q1 = a1/(nsites*2)
		q2 = a2/(nsites*2)
		Q = pd.DataFrame([q0, q1, q2]).T
		Q.columns = ['pop_0', 'pop_1', 'pop_2']
	elif n_anc == 4:
		a0 = arr[0::4] # should be views
		a1 = arr[1::4]
		a2 = arr[2::4]
		a3 = arr[3::4]
		q0 = a0/(nsites*2)
		q1 = a1/(nsites*2)
		q2 = a2/(nsites*2)
		q3 = a3/(nsites*2)
		Q = pd.DataFrame([q0, q1, q2, q3]).T
		Q.columns = ['pop_0', 'pop_1', 'pop_2', 'pop_3']

	return(Q)


def get_RMSD_Q(Q1, Q2):
	"""Return the RMSD between two sets of ancestry proportions."""
	assert(Q1.shape == Q2.shape)
	D = Q1-Q2
	SD = D*D
	MSD = SD.mean().mean()
	RMSD = np.sqrt(MSD)
	return(RMSD)


def max_la(vals, n_anc):
	"""
	Modifies the res array in place.
	Replaces each probabalistic la call with a categorical call
	for the ancestry with the highest posterior probability.
	Breaks ties by going to the lower numbered population
	"""
	idxs = np.arange(0, vals.shape[1]+n_anc, n_anc)
	for i in range(1, len(idxs)):
		b = vals[:, idxs[i-1]:idxs[i]]
		a = b.argmax(1, keepdims=True) # breaks ties by assigning the lower pop
		c = np.zeros_like(b)
		np.put_along_axis(c, a, 1, axis = 1)
		vals[:, idxs[i-1]:idxs[i]] = c
	return(vals)


def plot_ancestry_dosage(pred_dosage, start_index, n_anc, title, path=None, format='pdf', reference_dosage=None):
	"""    """
	colors = ['blue', 'orange', 'green', 'purple']

	fig, ax = plt.subplots(figsize = (12, n_anc*1.5), nrows=n_anc, sharex=True, sharey=True)
	f = []
	for i in range(n_anc):
		l, = ax[i].plot(pred_dosage[:, start_index+i], c=colors[i])
		f.append(l)
		if reference_dosage is not None:
			l, = ax[i].plot(reference_dosage[:, start_index+i], c=colors[i], alpha=.3, ls='--')

	plt.legend(f, [f'pop{p}' for p in range(n_anc)])



	fig.tight_layout()
	sns.despine(bottom=True)
	ax[0].set_title(title)
	ax[-1].set_xlabel('Site number ')
	if path:
		plt.savefig(path, dpi=300, format=format, bbox_inches='tight')
