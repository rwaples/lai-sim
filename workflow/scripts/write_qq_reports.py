import numpy as np
from common.utils import get_true_anc_dosage, get_ancestry_dosage, load_true_la, load_flare, load_mosaic, load_rfmix_fb, make_qq_report

# params
BCFTOOLS = str(snakemake.params.bcftools)
n_anc = int(snakemake.params.nsource)
# input
true_path = str(snakemake.input.true_la)
flare_path = str(snakemake.input.flare_la)
mosaic_path = str(snakemake.input.mosaic_la)
rfmix2_path = str(snakemake.input.rfmix2_la)

# output
flare_report = str(snakemake.output.flare_report)
mosaic_report = str(snakemake.output.mosaic_report)
rfmix_report = str(snakemake.output.rfmix_report)


true_anc_dosage = np.load(true_path)['arr_0']

flare_anc_dosage = np.load(flare_path)['arr_0']
flare_qq = make_qq_report(inferred_dosage=flare_anc_dosage, true_dosage=true_anc_dosage, nbins=400)
flare_qq.to_csv(flare_report, sep='\t', index=None, float_format='%.3f')
del flare_anc_dosage

try:
	mosaic_anc_dosage = np.load(mosaic_path)['arr_0']
	mosaic_qq = make_qq_report(inferred_dosage=mosaic_anc_dosage, true_dosage=true_anc_dosage, nbins=400)
	mosaic_qq.to_csv(mosaic_report, sep='\t', index=None, float_format='%.3f')
	del mosaic_anc_dosage
except ValueError:
	with open(mosaic_report, 'a'):  # Create file if does not exist
		pass

try:
	rfmix_anc_dosage = np.load(rfmix2_path)['arr_0']
	rfmix_qq = make_qq_report(inferred_dosage=rfmix_anc_dosage, true_dosage=true_anc_dosage, nbins=400)
	rfmix_qq.to_csv(rfmix_report, sep='\t', index=None, float_format='%.3f')
	del rfmix_anc_dosage
except ValueError:
	with open(rfmix_report, 'a'):  # Create file if does not exist
		pass
