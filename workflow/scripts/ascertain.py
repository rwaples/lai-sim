import tskit
import tszip
import numpy as np
from common.utils import strip_MAC, downsample_snps, get_allele_freqs, get_mean_allele_frequencies

ts_path = str(snakemake.input[0])
out_path = str(snakemake.output[0])
asc_MAC = int(snakemake.params.asc_MAC)
asc_MAF = float(snakemake.params.asc_MAF)
asc_maxsites = int(snakemake.params.asc_maxsites)
asc_seed = int(snakemake.params.asc_seed)
target_pop = int(snakemake.params.target_pop)

print(f'Start ascertinment: {out_path}')

ts = tszip.decompress(ts_path)

# apply MAC filter first
print(f'ascertainment MAC = {asc_MAC}')

ts = strip_MAC(ts, MAC=asc_MAC)

# apply MAF filter next
print(f'ascertainment MAF = {asc_MAF}')
if asc_MAF>0:
	af = get_allele_freqs(ts, pops=list(range(target_pop)))
	mean_af = get_mean_allele_frequencies(af)
	remove_index = np.where(mean_af<=asc_MAF)[0] # sites that failed maf filter
	to_remove = np.take(np.array([s.id for s in ts.sites()]), indices=remove_index)
	ts = ts.delete_sites(to_remove)


print(f'ascertainment max sites = {asc_maxsites}')
if asc_maxsites>0:
	if ts.num_sites > asc_maxsites:
		ts = downsample_snps(ts, nsnps=asc_maxsites, seed=asc_seed, fail=True)

tszip.compress(ts, out_path)
print(f'Done with ascertainment')
