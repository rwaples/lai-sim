mkdir ./programs


## BEAGEL 5.2
mkdir programs/BEAGLE
cd programs/RFmix
wget http://faculty.washington.edu/browning/beagle/beagle.29May21.d6d.jar


## RFMIX version 1 (1.5.4)
mkdir ./programs/RFmix
cd ./programs/RFmix
wget https://www.dropbox.com/s/cmq4saduh9gozi9/RFMix_v1.5.4.zip
unzip unzip ./RFMix_v1.5.4.zip
cd cd RFMix_v1.5.4
cd PopPhased
g++ -Wall -O3 -ftree-vectorize -fopenmp main.cpp getdata.cpp randomforest.cpp crfviterbi.cpp windowtosnp.cpp -o RFMix_PopPhased
cd ../TrioPhased/
g++ -Wall -O3 -ftree-vectorize -fopenmp main.cpp getdata.cpp randomforest.cpp crfviterbi.cpp windowtosnp.cpp -o RFMix_TrioPhased
cd ../../../..

## RFMIX version 2 (v2.03-r0)
mkdir ./programs/RFmix2
cd ./programs/RFmix2
git clone https://github.com/slowkoni/rfmix.git https://github.com/slowkoni/rfmix .
autoreconf --force --install # creates the configure script and all its dependencies
./configure                  # generates the Makefile
make
cd ../..


## ELAI
mkdir ./programs/elai
cd ./programs/elai
wget https://www.haplotype.org/download/elai-latest.tar.gz
tar -xvf ./elai-latest.tar.gz

## MOSAIC (https://maths.ucd.ie/~mst/MOSAIC/)
## requires R (4?)
mkdir ./programs/MOSAIC
cd ./programs/MOSAIC
conda activate ./r-env
wget https://csgitlab.ucd.ie/mst/mosaic/-/archive/master/mosaic-master.tar.gz
mkdir ./r-packages

#R CMD INSTALL -l ./r-packages mosaic-master.tar.gz

# from within R
install.packages("doParallel")
install.packages("Rcpp")
install.packages("ff")
install.packages("combinat")
install.packages("LaF")
install.packages("argparser")
install.packages("mosaic-master.tar.gz")

# done with R
tar -xvf mosaic-master.tar.gz
cd mosaic-master/

# MOSAIC
cd ./programs/MOSAIC/mosaic-master
Rscript mosaic.R Moroccan example_data/ -a 2 -n 2 -c 18:22
cd ../../..


# test installations
# notice change to "python2" different than manual
cd ./programs/RFmix/RFMix_v1.5.4
python2 RunRFMix.py TrioPhased ./TestData/alleles1.txt ./TestData/classes.txt ./TestData/markerLocationsChr1.txt -o outputTrioPhased
cd ../../..


## BEAGLE
java -jar programs/BEAGLE/beagle.29May21.d6d.jar \
gt=results/local_ancestry/AmericanAdmixture_4B11/AA_42/test_anal_1.genotypes.target_inds.vcf.gz \
ref=results/local_ancestry/AmericanAdmixture_4B11/AA_42/test_anal_1.genotypes.reference_inds.vcf.gz \
map=results/local_ancestry/AmericanAdmixture_4B11/AA_42/test_anal_1.genetic_map.txt \
out=~/temp
