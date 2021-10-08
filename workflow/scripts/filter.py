import tskit
import tszip
import numpy as np
from common.utils import strip_MAC, strip_adjacent_sites, downsample_snps

ts_path = str(snakemake.input[0])
ts_out = str(snakemake.output.tsz)
mapping_out = str(snakemake.output.node_mapping)

MAC_filter = int(snakemake.params.MAC_filter)
max_snps = int(snakemake.params.max_snps)
anal_seed = int(snakemake.params.anal_seed)

ts = tszip.decompress(ts_path)
simple_ts, node_mapping = ts.simplify(
	samples=np.where(ts.tables.nodes.asdict()['time']==0)[0],
	map_nodes=True,
	filter_populations=False,
	filter_individuals=True,
	filter_sites=True,
)
simple_ts = strip_MAC(simple_ts, MAC=MAC_filter)
simple_ts = strip_adjacent_sites(simple_ts, dist=1.5)
simple_ts = downsample_snps(simple_ts, nsnps=max_snps, seed=anal_seed, fail=True)

tszip.compress(simple_ts, ts_out)
np.save(file=mapping_out, arr=node_mapping)
