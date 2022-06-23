"""Load output from flare and save diploid ancestry dosages as compressed numpy array."""
import numpy as np
from common.utils import get_ancestry_dosage, load_flare
flare_path = str(snakemake.input.flare_la)
n_anc = int(snakemake.params.nsource)
sites_file = str(snakemake.params.sites_file)
out_path = str(snakemake.output)
BCFTOOLS = str(snakemake.params.BCFTOOLS)

flare_anc_dosage = get_ancestry_dosage(
	load_flare(
		flare_path,
		sites_file=sites_file,
		BCFTOOLS=BCFTOOLS
	),
	n_anc=n_anc
)
np.savez_compressed(out_path, flare_anc_dosage)

# can be loaded with:
# np.load(out_path)['arr_0']
