from common.utils import get_true_anc_dosage, get_ancestry_dosage, load_true_la, load_flare, load_mosaic, load_rfmix_fb, plot_ancestry_dosage

true_path = str(snakemake.input.true_la)
mosaic_path = str(snakemake.input.mosaic_la)
rfmix2_path = str(snakemake.input.rfmix2_la)
flare_path = str(snakemake.input.flare_la)
sites_file = str(snakemake.input.sites_file)


BCFTOOLS = str(snakemake.params.bcftools)
n_anc = int(snakemake.params.nsource)

path = str(snakemake.params.path)
format = str(snakemake.params.format)

nplot = 3  # number of diploids to plot

true_anc_dosage = get_true_anc_dosage(load_true_la(true_path), n_anc=n_anc)
rfmix_anc_dosage = get_ancestry_dosage(load_rfmix_fb(rfmix2_path), n_anc=n_anc)
mosaic_anc_dosage = get_ancestry_dosage(load_mosaic(mosaic_path), n_anc=n_anc)
flare_anc_dosage = get_ancestry_dosage(load_flare(flare_path, sites_file=sites_file, BCFTOOLS=BCFTOOLS), n_anc=n_anc)

for i in range(nplot):
	plot_ancestry_dosage(
		pred_dosage=rfmix_anc_dosage,
		start_index=i * n_anc,
		n_anc=n_anc,
		title=f'RFMix_v2 ancestry dosage \n ind {i}',
		path=path + f'.{i}.rfmix2.' + format,
		format=format,
		reference_dosage=true_anc_dosage
	)
	plot_ancestry_dosage(
		pred_dosage=mosaic_anc_dosage,
		start_index=i * n_anc,
		n_anc=n_anc,
		title=f'MOSAIC ancestry dosage \n ind {i}',
		path=path + f'.{i}.mosaic.' + format,
		format=format,
		reference_dosage=true_anc_dosage
	)
	plot_ancestry_dosage(
		pred_dosage=flare_anc_dosage,
		start_index=i * n_anc,
		n_anc=n_anc,
		title=f'flare ancestry dosage \n ind {i}',
		path=path + f'.{i}.flare.' + format,
		format=format,
		reference_dosage=true_anc_dosage
	)
