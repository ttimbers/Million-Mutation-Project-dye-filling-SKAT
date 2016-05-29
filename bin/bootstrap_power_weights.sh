# Shell script to iterate over for bootstrap power analysis
# Tiffany Timbers
# 2016-05-08

# export path to plink to $PATH
export PATH="/usr/bin/plink:$PATH"

# get container ID
cID=$(cat /proc/self/cgroup | grep "cpu:/" | sed 's/\([0-9]\):cpu:\/docker\///g')

## testing that we can make directories with the jobID
mkdir data/$cID

## Create list of randomly sampled strains (without replacement) & phenotype data from data/phenotype_amphid_dyf_dichotomous.csv
Rscript bin/create_random_samples.R data/phenotype_amphid_dyf_log_005.tsv \t TRUE 1 $1 data/$cID/temp_phenotype_amphid_dyf_log.tsv

## Create list of randomly selected strains from temp_phenotype_amphid_dyf.csv
awk '{print $1}' data/$cID/temp_phenotype_amphid_dyf_log.tsv | grep -h "^VC*" > data/$cID/temp_list_VCstrains_vcf.tsv

## Create a vcf file from these random selected strains, only variants from data/MMPfiltered.vcf
gunzip -c data/MMP.vcf.gz | perl bin/filter_MMP_variants.pl -input - -output data/$cID/temp_MMPcoding.vcf -strain data/$cID/temp_list_VCstrains_vcf.tsv -protein

## Create SSID file for SKAT analysis
Rscript bin/Make_SSID_file.R data/$cID/temp_MMPcoding.vcf data/$cID/temp_MMPcoding.SSID

## Create a filtered SSID file and vcf file for only variants from those genes which have
## a specified minimum number of alleles (we chose 7)
Rscript bin/create_reduced_variant_files.R data/$cID/temp_MMPcoding.vcf data/$cID/temp_MMPcoding.SSID 7 data/$cID/temp_MMPfiltered.vcf data/$cID/temp_MMPfiltered.SSID

## Create a custom weights for each variant based off of the type of mutation.
## Script Arguments:
## 1: name of .vcf file containing coding variants
## 2: weight to assign "knockout" mutations
## 3: weight to assign inframe deletions
## 4: weight to assign missense mutations
## 5: name of file to save variant weights to
Rscript bin/Make_custom_variant_weights.R data/$cID/temp_MMPfiltered.vcf 1 0.75 0.25 data/$cID/temp_MMP_SNP_WeightFile.txt

## Create binary plink files from the vcf file
#if [ ! -d "data/$cID" ]; then mkdir data/$cID; fi
plink --vcf data/$cID/temp_MMPfiltered.vcf --allow-extra-chr --no-fid --no-parents --no-sex --no-pheno --out data/$cID/temp_MMPfiltered

## Perform SKAT analysis
Rscript bin/do_SKAT_weights.R data/$cID/temp_MMPfiltered.fam data/$cID/temp_phenotype_amphid_dyf_log.tsv data/$cID data/$cID/temp_MMPfiltered.SSID data/$cID/temp_MMP_SNP_WeightFile.txt

## add N to dataframe
Rscript bin/add_column_to_tsv.R data/$cID/SKAT_weights.tsv $1 N data/$cID/SKAT_weights.tsv

## add container ID to dataframe
Rscript bin/add_column_to_tsv.R data/$cID/SKAT_weights.tsv $cID ID data/$cID/SKAT_weights.tsv

## clean up unneccesary files
rm data/$cID/temp* data/$cID/file*
