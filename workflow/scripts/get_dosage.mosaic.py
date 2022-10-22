"""Load output from mosaic and save diploid ancestry dosages as compressed numpy array."""
import numpy as np
from common.utils import get_ancestry_dosage, load_mosaic
mosaic_path = str(snakemake.input.mosaic_la)
n_anc = int(snakemake.params.nsource)
out_path = str(snakemake.output)

mosaic_anc_dosage = get_ancestry_dosage(
	load_mosaic(
		mosaic_path
	),
	n_anc=n_anc
)
np.savez_compressed(out_path, mosaic_anc_dosage)

# can be loaded with:
# np.load(out_path)['arr_0']
