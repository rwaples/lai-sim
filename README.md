# lai-sim
A Snakemake workflow for simulating admixture and evaluating local ancestry inference programs.

## Required external programs
* SLiM v3[link](https://messerlab.org/slim/)
* msprime [(link)](https://tskit.dev/msprime/docs/stable/intro.html)
* RFMix version 2 [(link)](https://github.com/slowkoni/rfmix)
* flare [(link)](https://github.com/browning-lab/flare)
* MOSAIC [(link)](https://maths.ucd.ie/~mst/MOSAIC/)
* BCFTOOLS [(link)](http://samtools.github.io/bcftools/howtos/index.html)
* BEAGLE [(link)](http://faculty.washington.edu/browning/beagle/beagle.html)
* add_err (no link provided, personal program of Brian L. Browning - browning@uw.edu)
* Snakemake [(link)](https://snakemake.readthedocs.io/en/stable/)
* conda [(link)](https://docs.conda.io/en/latest/)

### create and activate a conda environment

#### Snakemake suggest mamba (and it is much faster, but it breaks msprime on my machine)
#### So I am using conda below
#### uses mamba, rather than conda, as suggested by Snakemake
`mamba env create -f environment.yml --prefix ./env`
`conda env create -f environment.yml --prefix ./env`
#### activate conda environment (from base directory)
`conda activate ./env`

#### Install MOSAIC package
#### requires R v4
`mkdir ./programs/MOSAIC`
`cd ./programs/MOSAIC`
`wget https://maths.ucd.ie/~mst/MOSAIC/MOSAIC_1.3.9.tar.gz`
`# within R`
`install.packages("programs/MOSAIC/MOSAIC_1.3.9.tar.gz")`

### dry run of the test pipeline
#### defined in profiles/test/config.yaml
#### test analyses described in:
  - ./config/test.simulations.tsv
  - ./config/test.ascertainents.tsv,
  - ./config/test.analyses.tsv,
`snakemake --configfile ./profiles/test/config.yaml  --cores 11 --dry-run`

### run the test pipeline
`snakemake --configfile ./profiles/test/config.yaml --cores 11`

### results will be stored in:
  - ./results/reports/  (contains files with summaries of each analysis, for RMSD, R^2, QQ, etc.)
  - ./results/{model_name}/{simulation_name}/{ascertainment_name}/{analysis_name}/ (contains the simulated data used as input to LAI, as well as the LAI calls and evaluations vs truth)

