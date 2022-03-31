"""Load output from bmix and save diploid ancestry dosages as compressed numpy array."""
import numpy as np
from common.utils import get_ancestry_dosage, load_rfmix_fb
n_anc = int(snakemake.params.nsource)
rfmix2_path = str(snakemake.input.rfmix2_la)
out_path = str(snakemake.output)

rfmix2_anc_dosage = get_ancestry_dosage(
	load_rfmix_fb(
		rfmix2_path
	),
	n_anc=n_anc
)
np.savez_compressed(out_path, rfmix2_anc_dosage)

# can be loaded with:
# np.load(out_path)['arr_0']
