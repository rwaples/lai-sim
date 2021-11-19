import numpy as np
import pandas as pd
import pyreadr
import os

# params
BCFTOOLS = str(snakemake.params.bcftools)
n_anc = int(snakemake.params.nsource)
# input
true_path = str(snakemake.input.true_la)
mosaic_path = str(snakemake.input.mosaic_la)
rfmix2_path = str(snakemake.input.rfmix2_la)
bmix_path = str(snakemake.input.bmix_la)
# output
RMSD_path = str(snakemake.output.RMSD_path)
Q_true_path = str(snakemake.output.Q_true_path)
Q_bmix_path = str(snakemake.output.Q_bmix_path)
Q_mosaic_path = str(snakemake.output.Q_mosaic_path)
Q_rfmix_path = str(snakemake.output.Q_rfmix_path)

def get_ancestry_dosage(arr, n_anc):
    anc_dosage = np.zeros((arr.shape[0], int(arr.shape[1]/2)), dtype=np.half)
    if n_anc==3:
        assert (n_anc==3)
        a0 = arr[:, 0::3] # should be views
        a1 = arr[:, 1::3]
        a2 = arr[:, 2::3]
        anc_dosage[:, 0::3] = a0[:, ::2] + a0[:, 1::2]
        anc_dosage[:, 1::3] = a1[:, ::2] + a1[:, 1::2]
        anc_dosage[:, 2::3] = a2[:, ::2] + a2[:, 1::2]
    elif n_anc==4:
        assert (n_anc==4)
        a0 = arr[:, 0::4] # should be views
        a1 = arr[:, 1::4]
        a2 = arr[:, 2::4]
        a3 = arr[:, 3::4]
        anc_dosage[:, 0::4] = a0[:, ::2] + a0[:, 1::2]
        anc_dosage[:, 1::4] = a1[:, ::2] + a1[:, 1::2]
        anc_dosage[:, 2::4] = a2[:, ::2] + a2[:, 1::2]
        anc_dosage[:, 3::4] = a3[:, ::2] + a3[:, 1::2]
    return anc_dosage

def load_true_la(path):
	return np.load(path)['arr']


def get_true_anc_dosage(true_la, n_anc):
	hap1 = np.zeros((true_la.shape[0], int(true_la.shape[1]/2*n_anc)), dtype = 'int8')
	hap2 = np.zeros((true_la.shape[0], int(true_la.shape[1]/2*n_anc)), dtype = 'int8')
	aa = np.arange(true_la[:, ::2].shape[1])*n_anc+true_la[:, ::2]
	bb = np.arange(true_la[:, 1::2].shape[1])*n_anc+true_la[:, 1::2]
	np.put_along_axis(hap1, aa, 1, axis=1)
	np.put_along_axis(hap2, bb, 1, axis=1)
	return hap1+hap2


def r2_ancestry_dosage(true_dosage, pred_dosage, n_anc):
    per_anc = []
    for i in range(n_anc):
        per_anc.append(
            pearsonr(
                true_dosage[:,i::n_anc].ravel(),
                pred_dosage[:,i::n_anc].ravel()
            )[0]
        )
    per_ind = []
    for i in range(int(true_dosage.shape[1]/n_anc)):
        per_ind.append(
            pearsonr(
                true_dosage[:, i*n_anc:i*n_anc+n_anc].ravel(),
                pred_dosage[:, i*n_anc:i*n_anc+n_anc].ravel()
            )[0]
        )

    return(per_anc, per_ind)


## Load in the probablistic output of each method
def load_rfmix_fb(path):
	rfmix_res = pd.read_csv(path, sep='\t', comment='#')
	# expand out to each site
	rfmix_res = np.repeat(rfmix_res.iloc[:, 4:].values, [5], axis = 0)
	return rfmix_res

def load_bmix(path):
	csv_path = path.replace('.vcf.gz', '.csv')
	os.system(f"{BCFTOOLS} query -f '%CHROM, %POS, [%ANP1, %ANP2,]\\n' {path} > {csv_path}")
	bmix = pd.read_csv(csv_path, header=None)
	bmix = bmix.dropna(axis=1)
	return(bmix.iloc[:,2:].values)

def load_mosaic(path):
	print(path)
	mr = pyreadr.read_r(path)['arr'].astype(np.half)
	return(mr.to_numpy().T.reshape((mr.shape[2],-1), order='C'))


def get_Q(arr, n_anc):
    nsites = arr.shape[0]
    # avoid overflow and sum over sites
    arr = arr.astype(float).sum(0)
    if n_anc == 3:
        a0 = arr[0::3] # should be views
        a1 = arr[1::3]
        a2 = arr[2::3]
        q0 = (a0[0::2] + a0[1::2])/(nsites*4)
        q1 = (a1[0::2] + a1[1::2])/(nsites*4)
        q2 = (a2[0::2] + a2[1::2])/(nsites*4)
        Q = pd.DataFrame([q0, q1, q2]).T
        Q.columns = ['pop_0', 'pop_1', 'pop_2']
    elif n_anc == 4:
        a0 = arr[0::4] # should be views
        a1 = arr[1::4]
        a2 = arr[2::4]
        a3 = arr[3::4]
        q0 = (a0[0::2] + a0[1::2])/(nsites*4)
        q1 = (a1[0::2] + a1[1::2])/(nsites*4)
        q2 = (a2[0::2] + a2[1::2])/(nsites*4)
        q3 = (a3[0::2] + a3[1::2])/(nsites*4)
        Q = pd.DataFrame([q0, q1, q2, q3]).T
        Q.columns = ['pop_0', 'pop_1', 'pop_2', 'pop_3']

    return(Q)

def get_RMSD_Q(Q1, Q2):
    assert(Q1.shape == Q2.shape)
    D = Q1-Q2
    SD = D*D
    MSD = SD.mean().mean()
    RMSD = np.sqrt(MSD)
    return(RMSD)


true_anc_dosage = get_true_anc_dosage(load_true_la(true_path), n_anc=n_anc)
bmix_anc_dosage = get_ancestry_dosage(load_bmix(bmix_path), n_anc=n_anc)
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
    OUTFILE.write('\t'.join(['bmix', 'MOSAIC', 'RFMix2' ])  + '\n')
    OUTFILE.write('\t'.join([f'{x:0.4f}' for x in [rmsd_bmix, rmsd_mosaic, rmsd_rfmix]])  + '\n')

Q_true.to_csv(Q_true_path, index = None, sep = '\t', float_format='%0.4f')
Q_bmix.to_csv(Q_bmix_path, index = None, sep = '\t', float_format='%0.4f')
Q_mosaic.to_csv(Q_mosaic_path, index = None, sep = '\t', float_format='%0.4f')
Q_rfmix.to_csv(Q_rfmix_path, index = None, sep = '\t', float_format='%0.4f')
