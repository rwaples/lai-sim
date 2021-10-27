import numpy as np
import pandas as pd
import pyreadr
import os
from sklearn.metrics import r2_score

BCFTOOLS = str(snakemake.params.bcftools)
n_anc = int(snakemake.params.nsource)
true_path = str(snakemake.input.true_la)
mosaic_path = str(snakemake.input.mosaic_la)
rfmix2_path = str(snakemake.input.rfmix2_la)
bmix_path = str(snakemake.input.bmix_la)

R2_anc = str(snakemake.output.R2_anc)
R2_ind = str(snakemake.output.R2_ind)


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


def plot_ancestry_dosage(pred_dosage, start_index, reference_dosage=None):
	"""
	only works for 3 ancestries
	"""
	fig, ax = plt.subplots(figsize = (12, 4), nrows = 3, sharex=True, sharey=True)
	l0, = ax[0].plot(pred_dosage[:, start_index+0], c='b')
	l1, = ax[1].plot(pred_dosage[:, start_index+1], c='orange')
	l2, = ax[2].plot(pred_dosage[:, start_index+2], c='green')
	plt.legend([l0, l1, l2], ['pop0', 'pop1', 'pop2'])
	ax[0].set_title('Ancestry dosage')
	ax[2].set_xlabel('Site number ')

	if reference_dosage is not None:
		l0, = ax[0].plot(reference_dosage[:, start_index+0], c='b', alpha=.5, ls='--')
		l1, = ax[1].plot(reference_dosage[:, start_index+1], c='orange', alpha=.5, ls='--')
		l2, = ax[2].plot(reference_dosage[:, start_index+2], c='green', alpha=.5, ls='--')

	fig.tight_layout()
	sns.despine(bottom=True)



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
            r2_score(
                y_true=true_dosage[:,i::n_anc].ravel(),
                y_pred=pred_dosage[:,i::n_anc].ravel()
            )
        )
    per_ind = []
    for i in range(int(true_dosage.shape[1]/n_anc)):
        per_ind.append(
            r2_score(
                y_true=true_dosage[:, i*n_anc:i*n_anc+n_anc].ravel(),
                y_pred=pred_dosage[:, i*n_anc:i*n_anc+n_anc].ravel()
            )
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
	return mr.to_numpy().T.reshape((mr.shape[2],-1), order='C')



true_anc_dosage = get_true_anc_dosage(load_true_la(true_path), n_anc=n_anc)


rfmix_anc_dosage = get_ancestry_dosage(load_rfmix_fb(rfmix2_path), n_anc=n_anc)
rfmix_anc_r2, rfmix_ind_r2 = r2_ancestry_dosage(
	true_dosage=true_anc_dosage,
	pred_dosage=rfmix_anc_dosage,
	n_anc=n_anc
)


mosaic_anc_dosage = get_ancestry_dosage(load_mosaic(mosaic_path), n_anc=n_anc)
mosaic_anc_r2, mosaic_ind_r2 = r2_ancestry_dosage(
	true_dosage=true_anc_dosage,
	pred_dosage=mosaic_anc_dosage,
	n_anc=n_anc
)


bmix_anc_dosage = get_ancestry_dosage(load_bmix(bmix_path), n_anc=n_anc)
bmix_anc_r2, bmix_ind_r2 = r2_ancestry_dosage(
	true_dosage=true_anc_dosage,
	pred_dosage=bmix_anc_dosage,
	n_anc=n_anc
)


## Write R2 tables
with open(R2_anc, 'w') as OUTFILE:
	OUTFILE.write('\t'.join(['method'] + [f'anc_{x}' for x in range(n_anc)]) + '\n')
	OUTFILE.write('\t'.join(['rfmix2'] + [str(x) for x in rfmix_anc_r2])  + '\n')
	OUTFILE.write('\t'.join(['mosaic'] + [str(x) for x in mosaic_anc_r2])  + '\n')
	OUTFILE.write('\t'.join(['bmix'] + [str(x) for x in bmix_anc_r2])  + '\n')

with open(R2_ind, 'w') as OUTFILE:
	OUTFILE.write('\t'.join(['method'] + [f'ind_{x}' for x in range(len(bmix_ind_r2))]) + '\n')
	OUTFILE.write('\t'.join(['rfmix2'] + [str(x) for x in rfmix_ind_r2])  + '\n')
	OUTFILE.write('\t'.join(['mosaic'] + [str(x) for x in mosaic_ind_r2])  + '\n')
	OUTFILE.write('\t'.join(['bmix'] + [str(x) for x in bmix_ind_r2])  + '\n')
