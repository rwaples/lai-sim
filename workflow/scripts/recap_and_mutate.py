import tskit
import msprime
import stdpopsim
import tszip
import sys


ts_path = str(snakemake.input[0])
out_path = str(snakemake.output[0])
ancestral_Ne = int(snakemake.params.ancestral_Ne)
chr = str(snakemake.params.chr)
chr_len = float(snakemake.params.chr_len)
mutation_rate = float(snakemake.params.mutation_rate)
sim_seed = int(snakemake.params.sim_seed)

assert(mutation_rate < 1e-7) # make sure mutation rate is small

ts = tskit.load(ts_path)

species = stdpopsim.get_species("HomSap")
contig = species.get_contig(chr, length_multiplier=chr_len)

recapmap = msprime.RecombinationMap(
    positions=[0.0, ts.get_sequence_length()],
    rates=contig.recombination_map.get_rates(),
    num_loci=int(ts.get_sequence_length())
)

coalesced_ts = msprime.simulate(
    Ne=ancestral_Ne,
    from_ts=ts,
    random_seed=sim_seed,
    recombination_map=recapmap,
    population_configurations=[msprime.PopulationConfiguration() for p in ts.populations()],
    migration_matrix=None
)

mut_ts = msprime.mutate(tree_sequence=coalesced_ts, random_seed=sim_seed, rate=mutation_rate)
tszip.compress(mut_ts, out_path)
