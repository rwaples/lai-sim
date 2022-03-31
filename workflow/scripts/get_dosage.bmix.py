"""Load output from bmix and save diploid ancestry dosages as compressed numpy array."""
import numpy as np
from common.utils import get_ancestry_dosage, load_bmix
n_anc = int(snakemake.params.nsource)
bmix_path = str(snakemake.input.bmix_la)
sites_file = str(snakemake.input.sites_file)
out_path = str(snakemake.output)
BCFTOOLS = str(snakemake.params.bcftools)

bmix_anc_dosage = get_ancestry_dosage(
	load_bmix(
		bmix_path,
		sites_file=sites_file,
		BCFTOOLS=BCFTOOLS
	),
	n_anc=n_anc
)
np.savez_compressed(out_path, bmix_anc_dosage)

# can be loaded with:
# np.load(out_path)['arr_0']
