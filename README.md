# lai-sim
A Snakemake workflow for simulating admixture and evaluating local ancestry inference programs.

## Required external programs
* SLiM [link](https://messerlab.org/slim/)
* msprime [(link)](https://tskit.dev/msprime/docs/stable/intro.html)
* RFMix version 2 [(link)](https://github.com/slowkoni/rfmix)
* flare [(link)](https://github.com/browning-lab/flare)
* MOSAIC [(link)](https://maths.ucd.ie/~mst/MOSAIC/)
* BCFTOOLS [(link)](http://samtools.github.io/bcftools/howtos/index.html)
* BEAGLE [(link)](http://faculty.washington.edu/browning/beagle/beagle.html)
* add_err (no link provided)
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

### dry run of the test pipeline
`snakemake all -c1 --configfile profiles/personal/config.yaml --dry-run`

### run the full pipeline
`snakemake all -c1 --configfile profiles/cluster/config.yaml`

### view the directed acyclic graph (DAG) produced by Snakemake.
`snakemake all -c1 --configfile profiles/personal/config.yaml --rulegraph | dot | display`
