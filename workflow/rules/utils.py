import tskit


def strip_MAC(ts, MAC):
    """remove sites with a minor allele count <= MAC from the ts"""
    initial_sites = ts.num_sites
    num_samples = ts.num_samples
    tables = ts.dump_tables()
    tables.sites.clear()
    tables.mutations.clear()
    for tree in ts.trees():
        for site in tree.sites():
            assert len(site.mutations) == 1  # Only supports infinite sites muts.
            mut = site.mutations[0]
            if (tree.num_samples(mut.node) > MAC) and (tree.num_samples(mut.node) < num_samples-MAC):
                site_id = tables.sites.add_row(
                    position=site.position, ancestral_state=site.ancestral_state
                )
                tables.mutations.add_row(
                    site=site_id, node=mut.node, derived_state=mut.derived_state
                )
    final_sites = len(tables.sites)
    print(f"removed {initial_sites-final_sites} sites ({(initial_sites-final_sites)/(initial_sites):.0%}), {final_sites} sites remain")
    return(tables.tree_sequence())


def strip_adjacent(ts, dist = 1.5):
    """remove sites within dist bp apart from the ts"""

    initial_sites = ts.num_sites
    tables = ts.dump_tables()
    tables.sites.clear()
    tables.mutations.clear()

    prev_pos = 0
    for tree in ts.trees():
        for site in tree.sites():
            assert len(site.mutations) == 1  # Only supports infinite sites muts.
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
    print(f"removed {initial_sites-final_sites} sites ({(initial_sites-final_sites)/(initial_sites):.0%}), {final_sites} sites remain")
    return(tables.tree_sequence())


def downsample_snps(ts, nsnps, fail = False):
    """downsample a tree-sequence to a maximum number of snps
    if fail is True, will fail if there are not at least nsnps"""
    if ts.num_mutations == nsnps:
        return(ts)
    elif ts.num_mutations < nsnps:
        if fail:
            assert(False)
        else:
            return(ts)
    else:
        keep = frozenset(np.random.choice(a = ts.num_mutations, size = nsnps, replace = False))
        tables = ts.dump_tables()
        tables.sites.clear()
        tables.mutations.clear()
        mut_ix = 0
        for tree in ts.trees():
            for site in tree.sites():
                assert len(site.mutations) == 1  # Only supports infinite sites muts.
                if mut_ix in keep:
                    mut = site.mutations[0]
                    site_id = tables.sites.add_row(
                        position=site.position,
                        ancestral_state=site.ancestral_state)
                    tables.mutations.add_row(
                        site=site_id, node=mut.node, derived_state=mut.derived_state)
                mut_ix +=1
        return(tables.tree_sequence())


def get_local_ancestry(ts, admixture_time=20, nsample, per_rep = 12):

    l = [x for x in range(0, nsample, per_rep)]
    r = [x for x in range(per_rep, nsample+per_rep, per_rep)]
    dfs = []
    for i in range(len(l)):
        for pop in range(4):
            local = ts.tables.link_ancestors(
                # take the samples from each pop at time 0
                samples = np.intersect1d(
                    ts.samples(population = pop),
                    np.where(ts.tables.nodes.asdict()['time']==0)[0]
                )[l[i]:r[i]],
                ancestors = np.where(ts.tables.nodes.asdict()['time']==admixture_time)[0]
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
