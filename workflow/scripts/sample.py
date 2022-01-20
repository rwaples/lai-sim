import tskit
import tszip
import numpy as np
from numpy.random import default_rng
from common.utils import sample_inds

ts_path = str(snakemake.input[0])
out_path = str(snakemake.output[0])
nind_admixed = int(snakemake.params.nind_admixed)
nind_ref = str(snakemake.params.nind_ref)
admixture_time = int(snakemake.params.admixture_time)
anal_seed = int(snakemake.params.anal_seed)

nind_ref = np.array([int(x) for x in nind_ref.split(',')])


print(f'Start sampling: {out_path} (temp)')

ts = tszip.decompress(ts_path)

# select the samples
# assume the last popualtion is the admixed pop
pops = [p.id for p in ts.populations()]
ref_pops = pops[:-1]
admix_pops = pops[-1:]
to_take = np.array([], dtype = int)

# propagate the random seed here
rng = default_rng(anal_seed)
seeds = rng.bit_generator._seed_seq.spawn(len(pops))

i = 0
for pop in ref_pops:
	samples = sample_inds(ts, pop, nind_ref[i], seed = seeds[i])
	to_take = np.concatenate([to_take, samples])
	i+=1
for pop in admix_pops:
	samples = sample_inds(ts, pop, nind_admixed, seed = seeds[i])
	to_take = np.concatenate([to_take, samples])
del samples


# we also need to retain all the individuals present at admixture time
ancestors = np.where(ts.tables.nodes.asdict()['time']==admixture_time)[0]
to_take = np.concatenate([to_take, ancestors])
del ancestors
to_take.sort()


# simplyfy the TS
simple_ts = ts.simplify(
	samples=to_take,
	map_nodes=False,
	filter_populations=False,
	filter_individuals=False
)

# write the file back out
tszip.compress(simple_ts, out_path)

print(f'Done sampling')
