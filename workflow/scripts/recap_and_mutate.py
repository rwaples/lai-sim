import tskit
import msprime
import stdpopsim
import tszip
import sys


ts_path = str(snakemake.input[0])
out_path = str(snakemake.output[0])
ancestral_Ne = int(snakemake.params.ancestral_Ne)
chrom = str(snakemake.params.chrom)
contig_len = float(snakemake.params.contig_len)
mut_rate = float(snakemake.params.mut_rate)
seed = int(snakemake.params.seed)

assert(mut_rate < 1e-7) # make sure mutation rate is small

ts = tskit.load(ts_path)

species = stdpopsim.get_species("HomSap")
contig = species.get_contig(chrom, length_multiplier=contig_len)

recapmap = msprime.RecombinationMap(
    positions=[0.0, ts.get_sequence_length()],
    rates=contig.recombination_map.get_rates(),
    num_loci=int(ts.get_sequence_length())
)

coalesced_ts = msprime.simulate(
    Ne=ancestral_Ne,
    from_ts=ts,
    random_seed=seed,
    recombination_map=recapmap,
    population_configurations=[msprime.PopulationConfiguration() for p in ts.populations()],
    migration_matrix=None
)

mut_ts = msprime.mutate(tree_sequence = coalesced_ts, random_seed=seed, rate = mut_rate)
tszip.compress(mut_ts, out_path)
