os.system(f"{snakemake.config['PATHS']['BCFTOOLS']} query -f '%POS\\n' {os.path.join(base_path, 'genotypes.bcf')} > {os.path.join(base_path, 'site.positions')}")
positions = pd.read_csv(os.path.join(base_path, 'site.positions'), header=None)[0].values
