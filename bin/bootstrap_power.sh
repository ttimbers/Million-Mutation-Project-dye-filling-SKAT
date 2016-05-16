# Shell script to iterate over for bootstrap power analysis
# Tiffany Timbers
# 2016-05-08

# export path to plink to $PATH
export PATH="/usr/bin/plink:$PATH"

# get container ID
cID=$(cat /proc/self/cgroup | grep "cpu:/" | sed 's/\([0-9]\):cpu:\/docker\///g')

## testing that we can make directories with the jobID
mkdir data/$cID

## write a log file that matches input N with container ID
echo $cID','$1 > data/$cID/power.log

## Create list of randomly sampled strains (without replacement) & phenotype data from data/phenotype_amphid_dyf_dichotomous.csv
Rscript bin/create_random_samples.R data/phenotype_amphid_dyf_log.tsv \t TRUE 1 $1 data/$cID/temp_phenotype_amphid_dyf_log.tsv

## Create list of randomly selected strains from temp_phenotype_amphid_dyf.csv
awk '{print $1}' data/$cID/temp_phenotype_amphid_dyf_log.tsv | grep -h "^VC*" > data/$cID/temp_list_VCstrains_vcf.tsv

## Create a vcf file from these random selected strains, only variants from data/MMPfiltered.vcf
gunzip -c data/MMP.vcf.gz | perl bin/filter_MMP_variants.pl -input - -output data/$cID/temp_MMPcoding.vcf -strain data/$cID/temp_list_VCstrains_vcf.tsv -protein

## Create SSID file for SKAT analysis
Rscript bin/Make_SSID_file.R data/$cID/temp_MMPcoding.vcf data/$cID/temp_MMPcoding.SSID

## Create a filtered SSID file and vcf file for only variants from those genes which have
## a specified minimum number of alleles (we chose 7)
Rscript bin/create_reduced_variant_files.R data/$cID/temp_MMPcoding.vcf data/$cID/temp_MMPcoding.SSID 7 data/$cID/temp_MMPfiltered.vcf data/$cID/temp_MMPfiltered.SSID

## Create binary plink files from the vcf file
#if [ ! -d "data/$cID" ]; then mkdir data/$cID; fi
plink --vcf data/$cID/temp_MMPfiltered.vcf --allow-extra-chr --no-fid --no-parents --no-sex --no-pheno --out data/$cID/temp_MMPfiltered

## Perform SKAT analysis
Rscript bin/do_SKAT_no_weights.R data/$cID/temp_MMPfiltered.fam data/$cID/temp_phenotype_amphid_dyf_log.tsv data/$cID data/$cID/temp_MMPfiltered.SSID

## add N to dataframe
Rscript bin/add_column_to_tsv.R SKAT_pANDq_no_weights_results.txt $1 N SKAT_pANDq_no_weights_results.txt
