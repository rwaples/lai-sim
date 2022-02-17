# lai-sim
A Snakemake workflow for simulating admixture and evaluating local ancestry inference programs.

## Required external programs
	- SLiM [(link)](https://messerlab.org/slim/)
	- RFMix version 2 [(link)](https://github.com/slowkoni/rfmix)
	- bmix [(link)](https://github.com/browning-lab/bmix)
		BCFTOOLS [(link)](http://samtools.github.io/bcftools/howtos/index.html)
		BEAGLE [(link)](http://faculty.washington.edu/browning/beagle/beagle.html)
		add_err 
	- Snakemake and conda (see below)

# create the conda environment
## requires mamba, as suggested by Snakemake
`mamba env create -f environment.yml --prefix ./env`

# activate conda environment (from base directory)
`conda activate ./env`

# dry run the full pipeline
`snakemake all -c1 --configfile profiles/personal/config.yaml --dry-run`

# run the full pipeline with a single core
`snakemake all -c1 --configfile profiles/personal/config.yaml`

# view a graph of the directed acyclic graph (DAG) that is used for each analysis
`snakemake all -c1 --configfile profiles/personal/config.yaml --rulegraph | dot | display`

<!--- # patch MOSAIC to allow random number seeds
## needs to be run each time a new env is made
## patches the R code, allowing the passing of a random seed to the mosaic executable
`snakemake --cores 1 --force install_mosaic` --->

# to test the MOSAIC installation:
## run this cmd from programs/MOSAIC/MOSAIC
`Rscript ./mosaic.R simulated ./example_data/ -c 18:22 -n 3 -p "English Mandenka" --gens "30"`

# TODO
	- add seed to add-err.jar
