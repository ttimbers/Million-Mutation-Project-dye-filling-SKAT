# Shell script to iterate over for bootstrap power analysis
# Tiffany Timbers
# 2016-05-08

# export path to plink to $PATH
export PATH="/usr/bin/plink:$PATH"

# get container ID
cID=$(cat /proc/self/cgroup | grep "cpu:/" | sed 's/\([0-9]\):cpu:\/docker\///g')

## testing that the script runs
echo $PATH >> data/stdout.txt
echo $cID >> data/stdout.txt

## testing that we can make directories with the jobID
mkdir data/$cID

## testing that we can add .$PBS_JOBID to filenames
cp data/phenotype_amphid_dyf.csv data/phenotype_amphid_dyf.csv.$cID

## Create list of randomly sampled strains (without replacement) & phenotype data from data/phenotype_amphid_dyf_dichotomous.csv
Rscript bin/create_random_samples.R data/phenotype_amphid_dyf_dichotomous.csv \t TRUE 1 $1 data/temp_phenotype_amphid_dyf_dichotomous.csv.$cID

## Create list of randomly selected strains from temp_phenotype_amphid_dyf.csv
awk '{print $1}' data/temp_phenotype_amphid_dyf_dichotomous.csv.$cID | grep -h "^VC*" > data/temp_list_VCstrains_vcf.txt.$cID

## Create a vcf file from these random selected strains, only variants from data/MMPfiltered.vcf
gunzip -c data/MMP.vcf.gz | perl bin/filter_MMP_variants.pl -input - -output data/temp_MMPcoding.vcf.$cID -strain data/temp_list_VCstrains_vcf.txt.$cID -protein

## Create SSID file for SKAT analysis
Rscript bin/Make_SSID_file.R data/temp_MMPcoding.vcf.$cID data/temp_MMPcoding.SSID.$cID

## Create a filtered SSID file and vcf file for only variants from those genes which have
## a specified minimum number of alleles (we chose 7)
Rscript bin/create_reduced_variant_files.R data/temp_MMPcoding.vcf.$cID data/temp_MMPcoding.SSID.$cID 7 data/temp_MMPfiltered.vcf.$cID data/temp_MMPfiltered.SSID.$cID

## Create binary plink files from the vcf file
if [ ! -d "data/$cID" ]; then mkdir data/$cID; fi
plink --vcf data/temp_MMPfiltered.vcf.$cID --allow-extra-chr --no-fid --no-parents --no-sex --no-pheno --out data/power/temp_MMPfiltered

## Perform SKAT analysis
Rscript bin/do_SKAT_no_weights.R data/$cID/temp_MMPfiltered.fam data/temp_phenotype_amphid_dyf_dichotomous.csv.$cID data/$cID data/temp_MMPfiltered.SSID.$cID
