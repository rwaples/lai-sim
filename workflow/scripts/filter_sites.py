"""Select sites where local ancestry is recorded and
generate a simplified ts containing only these sites and the sampled nodes.
Notice this drop the ancestors retained to serve as local ancestry references. """
import tszip
import numpy as np
from common.utils import strip_MAC, strip_adjacent_sites, downsample_snps

ts_path = str(snakemake.input[0])
tsz_out = str(snakemake.output.tsz)
mapping_out = str(snakemake.output.node_mapping)

MAC_filter = int(snakemake.params.MAC_filter)
max_snps = int(snakemake.params.max_snps)
anal_seed = int(snakemake.params.anal_seed)


print(f'Start filter sites: {tsz_out}')

ts = tszip.decompress(ts_path)

simple_ts, node_mapping = ts.simplify(
	samples=np.where(ts.tables.nodes.asdict()['time'] == 0)[0],
	map_nodes=True,
	filter_populations=False,
	filter_individuals=True,
	filter_sites=True,
)
simple_ts = strip_MAC(simple_ts, MAC=MAC_filter)
simple_ts = strip_adjacent_sites(simple_ts, dist=1.5)
simple_ts = downsample_snps(simple_ts, nsnps=max_snps, seed=anal_seed, fail=False)

tszip.compress(simple_ts, tsz_out)
np.save(file=mapping_out, arr=node_mapping)

print('Done filter sites')
