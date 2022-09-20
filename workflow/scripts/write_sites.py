bcf = str(snakemake.input.bcf)
BCFTOOLS = str(snakemake.params.BCFTOOLS)
sites_path = str.output.sites

os.system(f"{BCFTOOLS} query -f '%POS\\n' {bcf} > {sites_path}")
