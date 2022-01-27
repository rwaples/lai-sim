import numpy as np
import pandas as pd
import os
from common.utils import max_la
from common.utils import get_true_anc_dosage, get_ancestry_dosage, load_true_la, load_bmix, load_mosaic, load_rfmix_fb, r2_ancestry_dosage


BCFTOOLS = str(snakemake.params.bcftools)
n_anc = int(snakemake.params.nsource)
true_path = str(snakemake.input.true_la)
mosaic_path = str(snakemake.input.mosaic_la)
rfmix2_path = str(snakemake.input.rfmix2_la)
bmix_path = str(snakemake.input.bmix_la)
sites_file = str(snakemake.input.sites_file)


R2_anc = str(snakemake.output.R2_anc)
R2_ind = str(snakemake.output.R2_ind)


true_anc_dosage = get_true_anc_dosage(load_true_la(true_path), n_anc=n_anc)


rfmix_anc_dosage = get_ancestry_dosage(load_rfmix_fb(rfmix2_path), n_anc=n_anc)
assert (len(rfmix_anc_dosage)-len(true_anc_dosage))<=5
rfmix_anc_r2, rfmix_ind_r2 = r2_ancestry_dosage(
	true_dosage=true_anc_dosage,
	# addressing possible uneven lengths due to RFMix2 only reporting every fifth site
	pred_dosage=rfmix_anc_dosage[:len(true_anc_dosage)],
	n_anc=n_anc
)
# get the r2 with max_like calls
rfmixML_anc_r2, rfmixML_ind_r2 = r2_ancestry_dosage(
	true_dosage=true_anc_dosage,
	# addressing possible uneven lengths due to RFMix2 only reporting every fifth site
	pred_dosage=max_la(rfmix_anc_dosage[:len(true_anc_dosage)], n_anc=n_anc),
	n_anc=n_anc
)
del rfmix_anc_dosage


mosaic_anc_dosage = get_ancestry_dosage(load_mosaic(mosaic_path), n_anc=n_anc)
mosaic_anc_r2, mosaic_ind_r2 = r2_ancestry_dosage(
	true_dosage=true_anc_dosage,
	pred_dosage=mosaic_anc_dosage,
	n_anc=n_anc
)


mosaicML_anc_r2, mosaicML_ind_r2 = r2_ancestry_dosage(
	true_dosage=true_anc_dosage,
	pred_dosage=max_la(mosaic_anc_dosage, n_anc=n_anc),
	n_anc=n_anc
)


del mosaic_anc_dosage


bmix_anc_dosage = get_ancestry_dosage(load_bmix(bmix_path, sites_file=sites_file, BCFTOOLS=BCFTOOLS), n_anc=n_anc)
bmix_anc_r2, bmix_ind_r2 = r2_ancestry_dosage(
	true_dosage=true_anc_dosage,
	pred_dosage=bmix_anc_dosage,
	n_anc=n_anc
)


bmixML_anc_r2, bmixML_ind_r2 = r2_ancestry_dosage(
	true_dosage=true_anc_dosage,
	pred_dosage=max_la(bmix_anc_dosage, n_anc=n_anc),
	n_anc=n_anc
)

del bmix_anc_dosage


## Write R2 tables
with open(R2_anc, 'w') as OUTFILE:
	OUTFILE.write('\t'.join(['method'] + [f'anc_{x}' for x in range(n_anc)]) + '\n')
	OUTFILE.write('\t'.join(['rfmix2'] + [f'{x:0.4f}' for x in rfmix_anc_r2])  + '\n')
	OUTFILE.write('\t'.join(['mosaic'] + [f'{x:0.4f}' for x in mosaic_anc_r2])  + '\n')
	OUTFILE.write('\t'.join(['bmix'] + [f'{x:0.4f}' for x in bmix_anc_r2])  + '\n')

	print(f'R^2 vs truth for LA calls:')
	print('\t'.join(['method'] + [f'anc_{x}' for x in range(n_anc)]))
	print('\t'.join(['rfmix2'] + [f'{x:0.4f}' for x in rfmix_anc_r2]))
	print('\t'.join(['mosaic'] + [f'{x:0.4f}' for x in mosaic_anc_r2]))
	print('\t'.join(['bmix'] + [f'{x:0.4f}' for x in bmix_anc_r2]))

with open(R2_ind, 'w') as OUTFILE:
	OUTFILE.write('\t'.join(['method'] + [f'ind_{x}' for x in range(len(bmix_ind_r2))]) + '\n')
	OUTFILE.write('\t'.join(['rfmix2'] + [f'{x:0.4f}' for x in rfmix_ind_r2])  + '\n')
	OUTFILE.write('\t'.join(['mosaic'] + [f'{x:0.4f}' for x in mosaic_ind_r2])  + '\n')
	OUTFILE.write('\t'.join(['bmix'] + [f'{x:0.4f}' for x in bmix_ind_r2])  + '\n')


with open(R2_anc+ ".ML", 'w') as OUTFILE:
	OUTFILE.write('\t'.join(['method'] + [f'anc_{x}' for x in range(n_anc)]) + '\n')
	OUTFILE.write('\t'.join(['rfmix2'] + [f'{x:0.4f}' for x in rfmixML_anc_r2])  + '\n')
	OUTFILE.write('\t'.join(['mosaic'] + [f'{x:0.4f}' for x in mosaicML_anc_r2])  + '\n')
	OUTFILE.write('\t'.join(['bmix'] + [f'{x:0.4f}' for x in bmixML_anc_r2])  + '\n')

	print(f'R^2 vs truth for [ML] LA calls:')
	print('\t'.join(['method'] + [f'anc_{x}' for x in range(n_anc)]))
	print('\t'.join(['rfmix2'] + [f'{x:0.4f}' for x in rfmixML_anc_r2]))
	print('\t'.join(['mosaic'] + [f'{x:0.4f}' for x in mosaicML_anc_r2]))
	print('\t'.join(['bmix'] + [f'{x:0.4f}' for x in bmixML_anc_r2]))

with open(R2_ind+ ".ML", 'w') as OUTFILE:
	OUTFILE.write('\t'.join(['method'] + [f'ind_{x}' for x in range(len(bmix_ind_r2))]) + '\n')
	OUTFILE.write('\t'.join(['rfmix2'] + [f'{x:0.4f}' for x in rfmixML_ind_r2])  + '\n')
	OUTFILE.write('\t'.join(['mosaic'] + [f'{x:0.4f}'for x in mosaicML_ind_r2])  + '\n')
	OUTFILE.write('\t'.join(['bmix'] + [f'{x:0.4f}' for x in bmixML_ind_r2])  + '\n')
