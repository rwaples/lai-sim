# lai-sim

# create the conda environment
## requires mamba, as suggested by Snakemake
`mamba env create -f environment.yml --prefix ./env`

# activate conda environment (from base directory)
`conda activate ./env`

# run the full pipeline with a single core
`snakemake --cores 1`

# view a graph of the DAG
`snakemake --dag | dot | display`

# re-install MOSAIC
snakemake --cores 1 --force install_mosaic

# test MOSAIC:
#### run from programs/MOSAIC/MOSAIC
Rscript ./mosaic.R simulated ./example_data/ -c 18:22 -n 3 -p "English Mandenka" --gens "30"


# TODO
	- ELAI
