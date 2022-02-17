from common.utils import get_true_anc_dosage, get_ancestry_dosage, load_true_la, load_bmix, load_mosaic, load_rfmix_fb, make_qq_report

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
bmix_report = str(snakemake.output.bmix_report)
mosaic_report = str(snakemake.output.mosaic_report)
rfmix_report = str(snakemake.output.rfmix_report)


true_anc_dosage = get_true_anc_dosage(load_true_la(true_path), n_anc=n_anc)
bmix_anc_dosage = get_ancestry_dosage(load_bmix(bmix_path, sites_file=sites_file, BCFTOOLS=BCFTOOLS), n_anc=n_anc)
mosaic_anc_dosage = get_ancestry_dosage(load_mosaic(mosaic_path), n_anc=n_anc)
rfmix_anc_dosage = get_ancestry_dosage(load_rfmix_fb(rfmix2_path), n_anc=n_anc)[:len(true_anc_dosage)]


bmix_qq = make_qq_report(inferred_dosage=bmix_anc_dosage, true_dosage=true_anc_dosage, nbins=13)
mosaic_qq = make_qq_report(inferred_dosage=mosaic_anc_dosage, true_dosage=true_anc_dosage, nbins=13)
rfmix_qq = make_qq_report(inferred_dosage=rfmix_anc_dosage, true_dosage=true_anc_dosage, nbins=13)

bmix_qq.to_csv(bmix_report, sep='\t', index=None, float_format='%.3f')
mosaic_qq.to_csv(mosaic_report, sep='\t', index=None, float_format='%.3f')
rfmix_qq.to_csv(rfmix_report, sep='\t', index=None, float_format='%.3f')
