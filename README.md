# lai-sim

# create the conda environment
## requires mamba, as suggested by Snakemake
`mamba env create -f environment.yml --prefix ./env`

# activate conda environment (from base directory)
`conda activate ./env`

# run the full pipeline with a single core
`snakemake --cores 1`

# view a graph of the directed acyclic graph (DAG) of the pipeline
`snakemake --dag | dot | display`

# patch MOSAIC to allow random number seeds
## needs to be run each time a new env is made
## patches the R code, allowing the passing of a random seed to the mosaic executable
`snakemake --cores 1 --force install_mosaic`

# test MOSAIC:
## run this cmd from programs/MOSAIC/MOSAIC
`Rscript ./mosaic.R simulated ./example_data/ -c 18:22 -n 3 -p "English Mandenka" --gens "30"`

# Benchmarking by snakemake benchmark:
`see: https://stackoverflow.com/questions/46813371/meaning-of-the-benchmark-variables-in-snakemake`


# TODO
	- remove indexing step from run_bmix (due to benchmarking)
	- generate a text table with running times for each analysis.
	- add seed to add-err.jar
	- fix issues with MOSAIC
