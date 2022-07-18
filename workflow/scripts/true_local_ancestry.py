import numpy as np
import pandas as pd
import tszip

from common.utils import get_local_ancestry_pop, get_la_mat_large

ancestry_ts = str(snakemake.input.ancestry_ts)
site_ts = str(snakemake.input.site_ts)
node_mapping = str(snakemake.input.node_mapping)
tracts_path = str(snakemake.output.tracts)
site_matrix_path = str(snakemake.output.site_matrix)
samples_path = str(snakemake.output.samples)
admixture_time = int(snakemake.params.admixture_time)
target_pop = int(snakemake.params.target_pop)
nind_admixed = int(snakemake.params.nind_admixed)


# local ancestry tracts
ts = tszip.decompress(ancestry_ts)
local_ancestry_df = get_local_ancestry_pop(
	ts,
	pop=target_pop,
	admixture_time=admixture_time,
	max_samples=nind_admixed * 2,
	per_rep=12
)
local_ancestry_df.to_hdf(
	tracts_path,
	key='local_ancestry',
	mode='w',
	complib='blosc:lz4',
	format='fixed'
)


# local ancestry per site
ts = tszip.decompress(site_ts)
mapping = np.load(file=node_mapping)
site_matrix, sample_order = get_la_mat_large(ts, df=local_ancestry_df, mapping=mapping)


np.savez_compressed(file=site_matrix_path, arr=site_matrix)
np.savetxt(fname=samples_path, X=sample_order, fmt='%i', delimiter="\t")
