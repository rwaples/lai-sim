"""calculates R2 score vs truth for each ancestry and each ind"""
import os
import numpy as np
from common.utils import max_la, r2_dosage_ancestry, r2_dosage_individual

n_anc = int(snakemake.params.nsource)
true_path = str(snakemake.input.true_dosage)
target_path = str(snakemake.input.target_dosage)

R2_anc = str(snakemake.output.R2_anc)
R2_ind = str(snakemake.output.R2_ind)

# load dosages
true_dosage = np.load(true_path)['arr_0'].astype(np.half)
# deal with proxy (empty) target dosage files.
try:
	target_dosage = np.load(target_path)['arr_0'].astype(np.half)
except ValueError:
	# check file is empty
	assert os.stat(target_path).st_size == 0
	# use an empty dosage as stand in.
	target_dosage = np.zeros_like(true_dosage, dtype=np.half)

assert (len(target_dosage) - len(true_dosage) <= 5)

anc_r2 = r2_dosage_ancestry(
	true_dosage=true_dosage,
	# addressing possible uneven lengths at the end by RFMix2
	pred_dosage=target_dosage[:len(true_dosage)],
	n_anc=n_anc
)

ind_r2 = r2_dosage_individual(
	true_dosage=true_dosage,
	# addressing possible uneven lengths at the end by RFMix2
	pred_dosage=target_dosage[:len(true_dosage)],
	n_anc=n_anc
)

# get the r2 with max_like calls
anc_r2ML = r2_dosage_ancestry(
	true_dosage=true_dosage,
	# addressing possible uneven lengths due to RFMix2 only reporting every fifth site
	pred_dosage=max_la(target_dosage[:len(true_dosage)], n_anc=n_anc),
	n_anc=n_anc
)

ind_r2ML = r2_dosage_individual(
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
