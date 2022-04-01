"""Load output from bmix and save diploid ancestry dosages as compressed numpy array."""
import numpy as np
from common.utils import get_true_anc_dosage, load_true_la
true_path = str(snakemake.input.true_la)
n_anc = int(snakemake.params.nsource)
out_path = str(snakemake.output)

true_anc_dosage = get_true_anc_dosage(
	load_true_la(true_path),
	n_anc=n_anc
)

np.savez_compressed(out_path, true_anc_dosage)

# can be loaded with:
# np.load(out_path)['arr_0']
