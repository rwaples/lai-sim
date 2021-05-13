import tskit
import tszip
from numpy.random import default_rng

MAC_filter = int(snakemake.params.MAC_filter)
max_snps = int(snakemake.params.max_snps)
random_seed = int(snakemake.params.random_seed)

random_seed = 42
rng = default_rng(random_seed)
