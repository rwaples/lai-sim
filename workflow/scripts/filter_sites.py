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
