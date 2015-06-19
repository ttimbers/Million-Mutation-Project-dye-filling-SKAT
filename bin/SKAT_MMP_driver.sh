## Assign dichotomous phenotype to strains. This script takes a filename of a .csv file 
##(which has 5 columns (Strain, Proportion of worms exhibiting phenotype 1, 
## Proportion of worms exhibiting phenotype 1, N, and Phenotype Summary),the name of the 
## control strain and the name of the first phenotype as well as the second phenotype to be 
## examined.
##
## Returns a .txt of 3 columns, named Strain, <phenotype1> and <phenotype2>. Entries for 
## <phenotype1> and <phenotype2> contain 1 if strain's phenotype diverges significantly 
## from wild-type (Fisher's exact test) and contain 1 if strain's phenotype diverges 
## significantly from wild-type, and 0 if it does not. A 5% FDR (Benjamini-Hochberg 
## procedure) is used to adjust for multiple comparisons.
Rscript bin/Assign_dichotomous_phenotype.R

## makes a file called phenotypes.txt from a file which lists strains in the first column 
## and phenotypes in all subsequent columns. Note input file must have a header listing
## the names of the phenotypes
## Script Arguments:
##    1. the name of the input file containing the phenotypes
##    2. the name of the output file
Rscript bin/Make_list_phenotype_file.R data/dyf_phenotypes_dichotomous.txt data/phenotypes.txt


## make list_VCstrains_vcf.txt file for strains assayed by extracting the first column 
## from a tab delimited text file to grab the strains from someone's tab delimited 
## phenotypes file:
awk '{print $1}' data/dyf_phenotypes_dichotomous.txt | grep -h "^VC*" > data/list_VCstrains_vcf.txt

## Create one merged .vcf file for only the strains that were assayed (listed in 
## a file called list_VCstrains_vcf.txt
## -protein argument creates a merged .vcf file with only variants which affect coding regions
gzcat data/MMP.vcf.gz | perl bin/filter_MMP_variants.pl -input - -output data/filteredMMP.vcf -strain data/list_VCstrains_vcf.txt -protein

## make a version of the merged .vcf file without the header
grep -v '^#' data/filteredMMP.vcf > data/MMP_non-syn_coding.vcf

## Make custom weights for each variant based off of the type of mutation. 
## Script Arguments:
## 1: name of .vcf file containing coding variants
## 2: weight to assign "knockout" mutations
## 3: weight to assign inframe deletions
## 4: weight to assign missense mutations
## 5: name of file to save variant weights to 
Rscript Make_custom_variant_weights.R data/MMP_non-syn_coding.vcf 1 0.75 0.25 MMP_SNP_WeightFile.txt

## Do SKAT analysis
Rscript do_SKAT.R data/MMP_non-syn_coding.vcf data/MMP_non-syn_coding_SSID.txt data/phenotypes.txt data/dyf_phenotypes_dichotomous.txt