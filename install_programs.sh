mkdir ./programs


# BCFTOOLS
cd programs
wget https://github.com/samtools/bcftools/releases/download/1.14/bcftools-1.14.tar.bz2
tar -xjvf bcftools-1.14.tar.bz2
cd bcftools-1.14/
./configure
make
cd ../..


## BEAGLE (v5.2 28Jun21)
mkdir programs/BEAGLE
cd programs/BEAGLE
wget https://faculty.washington.edu/browning/beagle/beagle.28Jun21.220.jar
cd ../..


# bmix
mkdir ./programs/bmix
cd ./programs/bmix
wget https://faculty.washington.edu/browning/bmix.jar
wget ../..


## RFMix2  (v2.03-r0)
mkdir programs/RFmix2
cd programs/RFmix2
git clone https://github.com/slowkoni/rfmix.git https://github.com/slowkoni/rfmix .
autoreconf --force --install # creates the configure script and all its dependencies
./configure                  # generates the Makefile
make
cd ../..


## MOSAIC (https://maths.ucd.ie/~mst/MOSAIC/)
## requires R v4
mkdir ./programs/MOSAIC
cd ./programs/MOSAIC
wget https://maths.ucd.ie/~mst/MOSAIC/MOSAIC_1.3.9.tar.gz
# within R
install.packages("programs/MOSAIC/MOSAIC_1.3.9.tar.gz")


<!--- # patch MOSAIC to allow random number seeds
## needs to be run each time a new env is made
## patches the R code, allowing the passing of a random seed to the mosaic executable
`snakemake --cores 1 --force install_mosaic` --->

### to test the MOSAIC installation:
run from programs/MOSAIC/MOSAIC
`Rscript ./mosaic.R simulated ./example_data/ -c 18:22 -n 3 -p "English Mandenka" --gens "30"`

### TODO
* add seed to add-err.jar




install.pacakge("reticulate")
# back in shell
git clone https://csgitlab.ucd.ie/mst/mosaic.git


# MOSAIC example
cd ./programs/MOSAIC/mosaic
Rscript mosaic.R Moroccan example_data/ -a 2 -n 2 -c 18:22
cd ../../..

## RFMIX version 1 (1.5.4)
mkdir programs/RFmix
cd programs/RFmix
wget https://www.dropbox.com/s/cmq4saduh9gozi9/RFMix_v1.5.4.zip
unzip ./RFMix_v1.5.4.zip
cd RFMix_v1.5.4
cd PopPhased
g++ -Wall -O3 -ftree-vectorize -fopenmp main.cpp getdata.cpp randomforest.cpp crfviterbi.cpp windowtosnp.cpp -o RFMix_PopPhased
cd ../TrioPhased/
g++ -Wall -O3 -ftree-vectorize -fopenmp main.cpp getdata.cpp randomforest.cpp crfviterbi.cpp windowtosnp.cpp -o RFMix_TrioPhased
cd ../../../..
# test installations
# notice change to "python2" different than manual
cd ./programs/RFmix/RFMix_v1.5.4
python2 RunRFMix.py TrioPhased ./TestData/alleles1.txt ./TestData/classes.txt ./TestData/markerLocationsChr1.txt -o outputTrioPhased
cd ../../..
## ELAI
mkdir ./programs/elai
cd ./programs/elai
wget https://www.haplotype.org/download/elai-latest.tar.gz
tar -xvf ./elai-latest.tar.gz
