def get_local_ancestry(ts, admixture_time, nsample, per_rep = 100):

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
