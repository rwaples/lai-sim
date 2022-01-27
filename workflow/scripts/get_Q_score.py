import numpy as np
import pandas as pd
import pyreadr
import os
from common.utils import get_true_anc_dosage, get_ancestry_dosage, load_true_la, load_bmix, load_mosaic, load_rfmix_fb, get_Q, get_RMSD_Q

# params
BCFTOOLS = str(snakemake.params.bcftools)
n_anc = int(snakemake.params.nsource)
# input
true_path = str(snakemake.input.true_la)
mosaic_path = str(snakemake.input.mosaic_la)
rfmix2_path = str(snakemake.input.rfmix2_la)
bmix_path = str(snakemake.input.bmix_la)
sites_file = str(snakemake.input.sites_file)
# output
RMSD_path = str(snakemake.output.RMSD_path)
Q_true_path = str(snakemake.output.Q_true_path)
Q_bmix_path = str(snakemake.output.Q_bmix_path)
Q_mosaic_path = str(snakemake.output.Q_mosaic_path)
Q_rfmix_path = str(snakemake.output.Q_rfmix_path)


true_anc_dosage = get_true_anc_dosage(load_true_la(true_path), n_anc=n_anc)
bmix_anc_dosage = get_ancestry_dosage(load_bmix(bmix_path, sites_file=sites_file, BCFTOOLS=BCFTOOLS), n_anc=n_anc)
mosaic_anc_dosage = get_ancestry_dosage(load_mosaic(mosaic_path), n_anc=n_anc)
rfmix_anc_dosage = get_ancestry_dosage(load_rfmix_fb(rfmix2_path), n_anc=n_anc)

Q_true = get_Q(true_anc_dosage, n_anc=n_anc)
Q_bmix = get_Q(bmix_anc_dosage, n_anc=n_anc)
Q_mosaic = get_Q(mosaic_anc_dosage, n_anc=n_anc)
Q_rfmix = get_Q(rfmix_anc_dosage, n_anc=n_anc)

rmsd_bmix = get_RMSD_Q(Q_bmix, Q_true)
rmsd_mosaic = get_RMSD_Q(Q_mosaic, Q_true)
rmsd_rfmix = get_RMSD_Q(Q_rfmix, Q_true)

## Write Q results tables
with open(RMSD_path, 'w') as OUTFILE:
	header = '\t'.join(['bmix', 'MOSAIC', 'RFMix2'])
	rmsd_line = '\t'.join([f'{x:0.4f}' for x in [rmsd_bmix, rmsd_mosaic, rmsd_rfmix]])
	OUTFILE.write(header + '\n')
	OUTFILE.write(rmsd_line  + '\n')
print(f'RMSD in Q values')
print(header)
print(rmsd_line)

Q_true.to_csv(Q_true_path, index = None, sep = '\t', float_format='%0.4f')
Q_bmix.to_csv(Q_bmix_path, index = None, sep = '\t', float_format='%0.4f')
Q_mosaic.to_csv(Q_mosaic_path, index = None, sep = '\t', float_format='%0.4f')
Q_rfmix.to_csv(Q_rfmix_path, index = None, sep = '\t', float_format='%0.4f')
