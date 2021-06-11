# lai-sim


# create the conda environment
## requires mamba
`mamba env create -f environment.yml --prefix ./env`

# re-install MOSAIC
snakemake --cores 1 --force install_mosaic

# test mosaic:
#### run from programs/MOSAIC/MOSAIC
Rscript ./mosaic.R simulated ./example_data/ -c 18:22 -n 3 -p "English Mandenka" --gens "30"


# TODO
	- figure out the 'distance' column in mosaic input
	-	summarize LAI output
	- ELAI 
