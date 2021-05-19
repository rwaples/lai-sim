import numpy as np
import pandas as pd
import collections


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
	return(samples)


def strip_MAC(ts, MAC):
    """remove sites with a minor allele count <= MAC from the ts"""
    initial_sites = ts.num_sites
    num_samples = ts.num_samples
    tables = ts.dump_tables()
    tables.sites.clear()
    tables.mutations.clear()
    for tree in ts.trees():
        for site in tree.sites():
            nmut = len(site.mutations)
            if(nmut==0):
                pass
            else:
                assert nmut== 1, f"{site}"  # Only supports infinite sites muts.
                mut = site.mutations[0]
                if (tree.num_samples(mut.node) > MAC) and (tree.num_samples(mut.node) < num_samples-MAC):
                    site_id = tables.sites.add_row(
                        position=site.position, ancestral_state=site.ancestral_state
                    )
                    tables.mutations.add_row(
                        site=site_id, node=mut.node, derived_state=mut.derived_state
                    )
    final_sites = len(tables.sites)
    print('MAC filter:')
    print(f"removed {initial_sites-final_sites} sites ({(initial_sites-final_sites)/(initial_sites):.0%}), {final_sites} sites remain")
    return(tables.tree_sequence())


def strip_adjacent(ts, dist):
    """remove sites within dist bp apart from the ts"""

    initial_sites = ts.num_sites
    tables = ts.dump_tables()
    tables.sites.clear()
    tables.mutations.clear()

    prev_pos = 0
    for tree in ts.trees():
        for site in tree.sites():
            nmut = len(site.mutations)
            if(nmut==0):
                pass
            else:
                assert nmut== 1, f"{site}"  # Only supports infinite sites muts.
                mut = site.mutations[0]
                pos = site.position
                if pos-prev_pos > dist:
                    site_id = tables.sites.add_row(
                        position=site.position, ancestral_state=site.ancestral_state
                    )
                    tables.mutations.add_row(
                        site=site_id, node=mut.node, derived_state=mut.derived_state
                    )
                prev_pos = pos
    final_sites = len(tables.sites)
    print('Adjacent site filter:')
    print(f"removed {initial_sites-final_sites} sites ({(initial_sites-final_sites)/(initial_sites):.0%}), {final_sites} sites remain")
    return(tables.tree_sequence())


def downsample_snps(ts, nsnps, seed, fail = False):
    """downsample a tree-sequence to a maximum number of snps
    if fail is True, will fail if there are not at least nsnps"""
    initial_sites = ts.num_sites
    if ts.num_mutations == nsnps:
        return(ts)
    elif ts.num_mutations < nsnps:
        if fail:
            assert False, "less than {nsnps} are present and fail is set to True"
        else:
            return(ts)
    else:
        rng = np.random.default_rng(seed)
        keep = frozenset(np.random.choice(a = ts.num_mutations, size = nsnps, replace = False))
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
                    assert nmut== 1, f"{site}"  # Only supports infinite sites muts.
                    if mut_ix in keep:
                        mut = site.mutations[0]
                        site_id = tables.sites.add_row(
                            position=site.position,
                            ancestral_state=site.ancestral_state)
                        tables.mutations.add_row(
                            site=site_id, node=mut.node, derived_state=mut.derived_state)
                    mut_ix +=1
        final_sites = len(tables.sites)
        print('Max SNPs filter:')
        print(f"removed {initial_sites-final_sites} sites ({(initial_sites-final_sites)/(initial_sites):.0%}), {final_sites} sites remain")
        return(tables.tree_sequence())


def get_local_ancestry(ts, admixture_time, per_batch):
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
    return(local_ancestry_df)



def get_la_mat(ts, df, mapping = None):
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

    return(la_mat, sample_order)


def make_ind_labels(export_ts):
	pop_of_sample = dict(zip(range(len(export_ts.tables.nodes)), export_ts.tables.nodes.population))
	nind_in_pop = collections.defaultdict(int)
	pops = [pop_of_sample[i] for i in export_ts.samples()]
	for p in pops:
		nind_in_pop[p] +=1
	for p in nind_in_pop:
		nind_in_pop[p] = int(nind_in_pop[p]/2)
	ind_labels = []
	for p in nind_in_pop:
		for i in range(1, nind_in_pop[p]+1):
			ind_labels.append(f'pop_{p}-ind_{i:04}')
	return(ind_labels)


def vcfheader(contig_id, contig_len, ind_labels):
    from datetime import date
    fileformat = '##fileformat=VCFv4.2\n'
    fileDate = f'##filedate={date.today().ctime()}\n'
    FILTER = '##FILTER=<ID=PASS,Description="All filters passed">\n'
    contig = f'##contig=<ID={contig_id},length={contig_len}>\n'
    FORMAT = '##FORMAT=<ID=LA,Number=1,Type=String,Description="Local ancestry">\n'
    header = fileformat+fileDate+FILTER+contig+FORMAT
    header = header +'\t'.join(['#CHROM','POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'])
    header = header + '\t' + '\t'.join(ind_labels) + '\n'
    return(header)
