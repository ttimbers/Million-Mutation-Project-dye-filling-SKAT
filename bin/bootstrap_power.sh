# Shell script to iterate over for bootstrap power analysis
# Tiffany Timbers
# 2016-05-08

## Create list of randomly sampled strains (without replacement) & phenotype data from data/phenotype_amphid_dyf_dichotomous.csv
Rscript bin/create_random_samples.R data/phenotype_amphid_dyf_dichotomous.csv \t TRUE 1 $1 data/temp_phenotype_amphid_dyf_dichotomous.csv.$PBS_JOBID

## Create list of randomly selected strains from temp_phenotype_amphid_dyf.csv
awk '{print $1}' data/temp_phenotype_amphid_dyf_dichotomous.csv.$PBS_JOBID | grep -h "^VC*" > data/temp_list_VCstrains_vcf.txt.$PBS_JOBID

## Create a vcf file from these random selected strains, only variants from data/MMPfiltered.vcf
gunzip -c data/MMP.vcf.gz | perl bin/filter_MMP_variants.pl -input - -output data/temp_MMPcoding.vcf.$PBS_JOBID -strain data/temp_list_VCstrains_vcf.txt.$PBS_JOBID -protein

## Create SSID file for SKAT analysis
Rscript bin/Make_SSID_file.R data/temp_MMPcoding.vcf.$PBS_JOBID data/temp_MMPcoding.SSID.$PBS_JOBID

## Create a filtered SSID file and vcf file for only variants from those genes which have
## a specified minimum number of alleles (we chose 7)
Rscript bin/create_reduced_variant_files.R data/temp_MMPcoding.vcf.$PBS_JOBID data/temp_MMPcoding.SSID.$PBS_JOBID 7 data/temp_MMPfiltered.vcf.$PBS_JOBID data/temp_MMPfiltered.SSID.$PBS_JOBID

## Create binary plink files from the vcf file
if [ ! -d "data/$PBS_JOBID" ]; then mkdir data/$PBS_JOBID; fi
plink --vcf data/temp_MMPfiltered.vcf.$PBS_JOBID --allow-extra-chr --no-fid --no-parents --no-sex --no-pheno --out data/power/temp_MMPfiltered

## Perform SKAT analysis
Rscript bin/do_SKAT_no_weights.R data/$PBS_JOBID/temp_MMPfiltered.fam data/temp_phenotype_amphid_dyf_dichotomous.csv.$PBS_JOBID data/$PBS_JOBID data/temp_MMPfiltered.SSID.$PBS_JOBID
