"""Load output from flare and save diploid ancestry dosages as compressed numpy array."""
import numpy as np
#import snakemake
from common.utils import get_ancestry_dosage, load_flare
flare_csv = str(snakemake.input.flare_csv)
flare_sites = str(snakemake.input.flare_sites)

n_anc = int(snakemake.params.nsource)
sites_file = str(snakemake.params.sites_file)
BCFTOOLS = str(snakemake.params.BCFTOOLS)

out_path = str(snakemake.output)

flare_anc_dosage = get_ancestry_dosage(
	load_flare(
		flare_csv,
		sites_file=sites_file,
		flare_sites=flare_sites,
		BCFTOOLS=BCFTOOLS
	),
	n_anc=n_anc
)
np.savez_compressed(out_path, flare_anc_dosage)

# can be loaded with:
# np.load(out_path)['arr_0']
