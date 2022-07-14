import numpy as np

flare_dosage = str(snakemake.input.dosage)
order_file = str(snakemake.input.order_file)
out_path = str(snakemake.output[0])
BCFTOOLS = str(snakemake.params.BCFTOOLS)
dosage = np.load(flare_dosage)['arr_0']


with open(order_file) as INFILE:
    tst = INFILE.readline().strip()
trt = tst.replace('##ANCESTRY=', '').strip('<>').split(',')
order = [x.split('=')[0] for x in trt]
check = [int(x.split('=')[1]) for x in trt]
assert check == list(range(len(check)))
n_anc = len(order)
n_ind = int(dosage.shape[1] / n_anc)
target_order = sorted(order)  # the order I want
reorder = np.array([target_order.index(x) for x in order])

dosage = dosage[:, reindex]
np.savez_compressed(out_path, dosage)
