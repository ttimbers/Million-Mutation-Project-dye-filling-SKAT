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
Rscript Assign_dichotomous_phenotype.R

## makes a file called phenotypes.txt from a file which lists strains in the first column 
## and phenotypes in all subsequent columns. Note input file must have a header listing
## the names of the phenotypes
## Script Arguments:
##    1. the name of the input file containing the phenotypes
##    2. the name of the output file
Rscript Make_list_phenotype_file.R dyf_phenotypes_dichotomous.txt phenotypes.txt


## make list_VCstrains_vcf.txt file for strains assayed by extracting the first column 
## from a tab delimited text file to grab the strains from someone's tab delimited 
## phenotypes file:
awk '{print $1}' dyf_phenotypes_dichotomous.txt | grep -h "^VC*" > vcf_files/list_VCstrains_vcf.txt

## Create one merged .vcf file for all strains that were assayed (listed in 
## a file called list_VCstrains_vcf.txt
##
## note: this .pl script also requires the use of a text file called 
## gknames.txt which contains all variant names of the MMP library, as 
## well as the VC*.combined.vcf files to combine into the merged file.
## These all should be in the same directory
cd vcf_files
./create_gwas_vcf_MMP.pl > ../MMP.vcf
cd ..

## create a merged .vcf file with only variants which affect coding regions
grep -v "GRANT=0" MMP.vcf > MMP_non-syn.vcf
grep -h "^#" MMP_non-syn.vcf > MMP_non-syn_coding.vcf
grep -h "coding" MMP_non-syn.vcf >> MMP_non-syn_coding_temp.vcf
grep -h -v "Noncoding" MMP_non-syn_coding_temp.vcf >> MMP_non-syn_coding.vcf
grep -h "GF=intron_splicing" MMP_non-syn.vcf >> MMP_non-syn_coding.vcf
rm MMP.vcf
rm MMP_non-syn.vcf

## make a version of the merged .vcf file without the header
grep -v '^#' MMP_non-syn_coding.vcf > MMP_non-syn_coding.vcf

## Make custom weights for each variant based off of the type of mutation. 
## Script Arguments:
## 1: name of .vcf file containing coding variants
## 2: weight to assign "knockout" mutations
## 3: weight to assign inframe deletions
## 4: weight to assign missense mutations
## 5: name of file to save variant weights to 
Rscript Make_custom_variant_weights.R MMP_non-syn_coding.vcf 1 0.75 0.25 MMP_SNP_WeightFile.txt

## Do SKAT analysis
## Rscript do_SKAT.R MMP_non-syn_coding.vcf MMP_non-syn_coding_SSID.txt phenotypes.txt dyf_phenotypes_dichotomous.txt

## Make snp set ID file (SSID) to run SKAT analysis
Rscript Make_SSID_file.R

## Make plink files from merged .vcf file to run SKAT analysis
Rscript Make_plink_files.R

## Do SKAT analysis
Rscript SKAT_logistic_strain.R