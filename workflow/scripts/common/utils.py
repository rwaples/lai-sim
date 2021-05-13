def sample_inds(ts, pop_id, nind, seed):
	"""return the haploid sample ids representing sampling nind individuals from pop_id in ts"""
	rng = np.random.default_rng(seed)
	hap_samples = ts.samples(population = pop_id)
	# sample from the first haploids of each ind
	take = rng.choice(hap_samples[::2], nind, replace=False)
	samples = np.empty(nind*2, dtype=int)
	samples[0::2] = take
	samples[1::2] = take+1
	return(samples)
