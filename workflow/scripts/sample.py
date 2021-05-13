import tskit
import tszip
import numpy as np
from numpy.random import default_rng
from common.utils import sample_inds

ts_path = str(snakemake.input[0])
out_path = str(snakemake.output[0])
nind_admixed = int(snakemake.params.nind_admixed)
nind_ref = int(snakemake.params.nind_ref)
random_seed = int(snakemake.params.random_seed)

rng = default_rng(random_seed)
# figure out a way to propagate the random seeds here


ts = tszip.decompress(ts_path)

# select the samples
# assume the last popualtions is the admixed pop
pops = [p.id for p in ts.populations()]
ref_pops = pops[:-1]
admix_pops = pops[-1:]
to_take = np.array([], dtype = int)
for pop in ref_pops:
	samples = sample_inds(ts, pop, nind_ref, seed = )
	to_take = np.concatenate([to_take, samples])
for pop in admix_pops:
	samples = sample_inds(ts, pop, nind_admixed)
	to_take = np.concatenate([to_take, samples])
to_take.sort()


# simplyfy the TS
simp_ts = ts.simplify(
	samples=to_take,
	map_nodes=False,
	filter_populations=False
)

# write the file back out
tszip.compress(out_path)

sample_inds
