import numpy as np
import pandas as pd
import pyreadr
import os.path
from common.utils import get_Q, get_RMSD_Q

n_anc = int(snakemake.params.nsource)
true_path = str(snakemake.input.true_dosage)
target_path = str(snakemake.input.target_path)
# output
Q_true_path = str(snakemake.output.Q_true)
Q_path = str(snakemake.output.Q)
RMSD_path = str(snakemake.output.RMSD)

true_dosage = np.load(true_path)['arr_0']
target_dosage = np.load(target_path)['arr_0']

Q_true = get_Q(true_dosage, n_anc=n_anc)
Q_target = get_Q(target_dosage, n_anc=n_anc)
RMSD = get_RMSD_Q(Q_target, Q_true)

Q_target.to_csv(Q_path, index=None, sep='\t', float_format='%0.4f')
if os.path.exists(Q_true_path):
	pass
else:
	Q_true.to_csv(Q_true_path, index=None, sep='\t', float_format='%0.4f')

# Write RMSD
with open(RMSD_path, 'w') as OUTFILE:
	OUTFILE.write(f'{RMSD:0.4f}' + '\n')
