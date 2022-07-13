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


# FLARE
mkdir ./programs/flare
cd ./programs/flare
wget https://faculty.washington.edu/browning/bmix.jar  # TODO: needs to be updated
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

### to test the MOSAIC installation:
run from programs/MOSAIC/MOSAIC
`Rscript ./mosaic.R simulated ./example_data/ -c 18:22 -n 3 -p "English Mandenka" --gens "30"`

