import tskit
import msprime
#import stdpopsim
import tszip
import numpy as np
from common.utils import strip_MAC


ts_path = str(snakemake.input.ts_path)
out_path = str(snakemake.output[0])
ancestral_Ne = int(snakemake.params.ancestral_Ne)
# chr = str(snakemake.params.chr)
# chr_len = float(snakemake.params.chr_len)
rec_map_path = str(snakemake.params.rec_map_path)
#max_bp = int(snakemake.params.max_bp)

mutation_rate = float(snakemake.params.mutation_rate)
sim_seed = int(snakemake.params.sim_seed)
admixture_time = int(snakemake.params.admixture_time)
assert(mutation_rate <= 1e-7)  # make sure mutation rate is not crazy high

ts = tskit.load(ts_path)


rmap = msprime.RateMap.read_hapmap(rec_map_path)
#recapmap = rmap.slice(0, max_bp + 1, trim=True)
recapmap = rmap

print('ts_length', ts.get_sequence_length())
#print('recapmap', recapmap.right[-1])
print('rmap', rmap.right[-1])


# species = stdpopsim.get_species("HomSap")
# contig = species.get_contig(chr, length_multiplier=chr_len)


# recombination map used for recapitation
# recapmap = msprime.RecombinationMap(
# positions=[0.0, ts.get_sequence_length()],
# rates=contig.recombination_map.get_rates(),
# num_loci=int(ts.get_sequence_length())
# )


# recapitate
# uses coalescent
coalesced_ts = msprime.simulate(
	Ne=ancestral_Ne,
	from_ts=ts,
	random_seed=sim_seed,
	recombination_map=recapmap,
	population_configurations=[msprime.PopulationConfiguration() for p in ts.populations()],
	migration_matrix=None
)

# Edit the ts so that only nodes at time=0 are marked as 'samples'
tables = coalesced_ts.tables
newnodes = tables.nodes.copy()

flags = tables.nodes.flags
assert flags.sum() == coalesced_ts.num_samples

flags[tables.nodes.time == admixture_time] = 0
assert flags.sum() == np.intersect1d(
	coalesced_ts.samples(),
	np.where(tables.nodes.asdict()['time'] == 0)[0]
).size

newnodes.set_columns(
	time=newnodes.time,
	flags=flags,
	population=newnodes.population,
	individual=newnodes.individual,
	metadata=newnodes.metadata,
	metadata_offset=newnodes.metadata_offset,
)

tables.nodes.clear()

for row in newnodes:
	tables.nodes.add_row(
		flags=row.flags,
		time=row.time,
		population=row.population,
		individual=row.individual,
		metadata=row.metadata
	)

coalesced_ts = tables.tree_sequence()
del newnodes
del tables
del flags


# add mutations
mut_ts = msprime.mutate(tree_sequence=coalesced_ts, random_seed=sim_seed, rate=mutation_rate)
# remove monomorphic sites
# monomorphic sites are present in the ts because some non-ancestors present in the ts.
mut_ts = strip_MAC(mut_ts, MAC=0)
# compress the ts with tszip
tszip.compress(mut_ts, out_path)
