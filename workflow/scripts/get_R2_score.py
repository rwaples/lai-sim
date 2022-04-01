"""calculates R2 score vs truth for each ancestry and each ind"""
import numpy as np
from common.utils import max_la, r2_ancestry_dosage

n_anc = int(snakemake.params.nsource)
true_path = str(snakemake.input.true_dosage)
target_path = str(snakemake.input.target_dosage)

R2_anc = str(snakemake.output.R2_anc)
R2_ind = str(snakemake.output.R2_ind)

# load dosages
true_dosage = np.load(true_path)['arr_0']
target_dosage = np.load(target_path)['arr_0']
assert (len(target_dosage) - len(true_dosage) <= 5)

anc_r2, ind_r2 = r2_ancestry_dosage(
	true_dosage=true_dosage,
	# addressing possible uneven lengths at the end by RFMix2
	pred_dosage=target_dosage[:len(true_dosage)],
	n_anc=n_anc
)

# get the r2 with max_like calls
anc_r2ML, ind_r2ML = r2_ancestry_dosage(
	true_dosage=true_dosage,
	# addressing possible uneven lengths due to RFMix2 only reporting every fifth site
	pred_dosage=max_la(target_dosage[:len(true_dosage)], n_anc=n_anc),
	n_anc=n_anc
)

# Write R2 tables
with open(R2_anc, 'w') as OUTFILE:
	OUTFILE.write('\t'.join([f'anc_{x}' for x in range(n_anc)]) + '\n')
	OUTFILE.write('\t'.join([f'{x:0.4f}' for x in anc_r2]) + '\n')

with open(R2_ind, 'w') as OUTFILE:
	OUTFILE.write('\t'.join([f'ind_{x}' for x in range(len(ind_r2))]) + '\n')
	OUTFILE.write('\t'.join([f'{x:0.4f}' for x in ind_r2]) + '\n')


with open(R2_anc + ".ML", 'w') as OUTFILE:
	OUTFILE.write('\t'.join([f'anc_{x}' for x in range(n_anc)]) + '\n')
	OUTFILE.write('\t'.join([f'{x:0.4f}' for x in anc_r2ML]) + '\n')

with open(R2_ind + ".ML", 'w') as OUTFILE:
	OUTFILE.write('\t'.join([f'ind_{x}' for x in range(len(ind_r2ML))]) + '\n')
	OUTFILE.write('\t'.join([f'{x:0.4f}' for x in ind_r2ML]) + '\n')
