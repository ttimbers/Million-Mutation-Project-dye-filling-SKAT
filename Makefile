## Tiffany A. Timbers
## July 25, 2015
##
## Makefile pipeline for SKAT analysis of genomic and phenotype data of Million Mutation 
## Project strains dye-filling phenotypes reported in Timbers et al.
##
## Dependencies: Bash Shell, Make, Perl, plink v1.90b1g (assumes that you have made plink 
## executable and put it in your Bash Shell's $PATH), R and R packages 
## SKAT (version 0.95), stringr, fdrtool


## Final files to be generated from SKAT analysis
all: data/amphid_dyf/SKAT_no_weights_results.txt data/amphid_dyf/SKAT_weights_results.txt data/amphid_dyf/SKAT_pANDq_no_weights_results.txt data/amphid_dyf/SKAT_pANDq_weights_results.txt data/phasmid_dyf/SKAT_no_weights_results.txt data/phasmid_dyf/SKAT_weights_results.txt data/phasmid_dyf/SKAT_pANDq_no_weights_results.txt data/phasmid_dyf/SKAT_pANDq_weights_results.txt

##======================================================================================
## Analysis of amphid dye-filling phenotype
##======================================================================================

## Assign dichotomous amphid dye-filling phenotype to strains. This script takes 2 
## arguments: (1) a comma delimited .csv file (which has 5 columns: Strain, 
## phenotype, N), and (2) output filename.
##  
## Returns a tab delimited .csv file (which has 4 columns: strain, phenotype, pvals, 
## fdr.bh). The phenotype column contains 1 if strain's phenotype diverges significantly 
## from wild-type (Fisher's exact test) and contain 1 if strain's phenotype diverges 
## significantly from wild-type, and 0 if it does not. A 5% FDR (Benjamini-Hochberg 
## procedure) is used to adjust for multiple comparisons.
data/phenotype_amphid_dyf_dichotomous.csv: bin/Assign_dichotomous_phenotype.R data/phenotype_amphid_dyf.csv
	Rscript bin/Assign_dichotomous_phenotype.R data/phenotype_amphid_dyf.csv data/phenotype_amphid_dyf_dichotomous.csv

## Extracts the strain column from the phenotype_amphid_dyf_dichotomous.csv file created 
## above and saves it as list_VCstrains_vcf.txt
data/list_VCstrains_vcf.txt: data/phenotype_amphid_dyf_dichotomous.csv
	awk '{print $$1}' data/phenotype_amphid_dyf_dichotomous.csv | grep -h "^VC*" > data/list_VCstrains_vcf.txt

## Creates one merged .vcf file for only the strains that were assayed (listed in a file
## called list_VCstrains_vcf.txt. Uses -protein argument creates a merged .vcf file with
## only variants which affect coding regions	
data/MMPcoding.vcf: data/MMP.vcf.gz bin/filter_MMP_variants.pl data/list_VCstrains_vcf.txt	
	gzcat data/MMP.vcf.gz | perl bin/filter_MMP_variants.pl -input - -output data/MMPcoding.vcf -strain data/list_VCstrains_vcf.txt -protein

## Create a version of the merged .vcf file without the header	
## data/MMPcoding_noHeader.vcf: data/MMPcoding.vcf
##	grep -v '^#' data/MMPcoding.vcf > data/MMPcoding_noHeader.vcf

## Create SSID file for SKAT analysis
data/MMPcoding.SSID: bin/Make_SSID_file.R data/MMPcoding.vcf
	Rscript bin/Make_SSID_file.R data/MMPcoding.vcf data/MMPcoding.SSID

## Create a filtered SSID file and vcf file for only variants from those genes which have
## a specified minimum number of alleles (we chose 6)
data/MMPfiltered.vcf data/MMPfiltered.SSID: bin/create_reduced_variant_files.R data/MMPcoding.vcf data/MMPcoding.SSID
	rscript bin/create_reduced_variant_files.R data/MMPcoding.vcf data/MMPcoding.SSID 6 data/MMPfiltered.vcf data/MMPfiltered.SSID

## Create a custom weights for each variant based off of the type of mutation. 
## Script Arguments:
## 1: name of .vcf file containing coding variants
## 2: weight to assign "knockout" mutations
## 3: weight to assign inframe deletions
## 4: weight to assign missense mutations
## 5: name of file to save variant weights to 	
data/MMP_SNP_WeightFile.txt: bin/Make_custom_variant_weights.R data/MMPfiltered.vcf
	Rscript bin/Make_custom_variant_weights.R data/MMPfiltered.vcf 1 0.75 0.25 data/MMP_SNP_WeightFile.txt	

## Create binary plink files for amphid phenotype from filtered .vcf file
data/amphid_dyf data/amphid_dyf/MMPfiltered.fam data/amphid_dyf/MMPfiltered.bim data/amphid_dyf/MMPfiltered.bed data/amphid_dyf/MMPfiltered.log: data/MMPfiltered.vcf
	mkdir data/amphid_dyf
	plink --vcf data/MMPfiltered.vcf --allow-extra-chr --no-fid --no-parents --no-sex --no-pheno --out data/amphid_dyf/MMPfiltered

## Perform SKAT analysis
data/amphid_dyf/SKAT_no_weights_results.txt data/amphid_dyf/SKAT_weights_results.txt data/amphid_dyf/SKAT_pANDq_no_weights_results.txt data/amphid_dyf/SKAT_pANDq_weights_results.txt: bin/do_SKAT.R data/amphid_dyf/MMPfiltered.fam data/MMP_SNP_WeightFile.txt data/MMPfiltered.SSID
	Rscript bin/do_SKAT.R data/amphid_dyf/MMPfiltered.fam data/phenotype_amphid_dyf_dichotomous.csv data/amphid_dyf data/MMPfiltered.SSID data/MMP_SNP_WeightFile.txt


##======================================================================================
## Analysis of phasmid dye-filling phenotype
##======================================================================================

## Assign dichotomous phasmid dye-filling phenotype to strains. This script takes 2 
## arguments: (1) a comma delimited .csv file (which has 5 columns: Strain, 
## phenotype, N), and (2) output filename.
##  
## Returns a tab delimited .csv file (which has 4 columns: strain, phenotype, pvals, 
## fdr.bh). The phenotype column contains 1 if strain's phenotype diverges significantly 
## from wild-type (Fisher's exact test) and contain 1 if strain's phenotype diverges 
## significantly from wild-type, and 0 if it does not. A 5% FDR (Benjamini-Hochberg 
## procedure) is used to adjust for multiple comparisons.
data/phenotype_phasmid_dyf_dichotomous.csv: bin/Assign_dichotomous_phenotype.R data/phenotype_phasmid_dyf.csv
	Rscript bin/Assign_dichotomous_phenotype.R data/phenotype_phasmid_dyf.csv data/phenotype_phasmid_dyf_dichotomous.csv

## Create binary plink files for phasmid phenotype from filtered .vcf file
data/phasmid_dyf data/phasmid_dyf/MMPfiltered.fam data/phasmid_dyf/MMPfiltered.bim data/phasmid_dyf/MMPfiltered.bed data/phasmid_dyf/MMPfiltered.log: data/MMPfiltered.vcf
	mkdir data/phasmid_dyf
	plink --vcf data/MMPfiltered.vcf --allow-extra-chr --out data/phasmid_dyf/MMPfiltered

## Perform SKAT analysis (use SSID file made above)
data/phasmid_dyf/SKAT_no_weights_results.txt data/phasmid_dyf/SKAT_weights_results.txt data/phasmid_dyf/SKAT_pANDq_no_weights_results.txt data/phasmid_dyf/SKAT_pANDq_weights_results.txt: bin/do_SKAT.R data/phasmid_dyf/filteredMMP.fam data/MMP_SNP_WeightFile.txt data/MMPfiltered.SSID data/phenotype_phasmid_dyf_dichotomous.csv
	Rscript bin/do_SKAT.R data/phasmid_dyf/MMPfiltered.fam data/phenotype_phasmid_dyf_dichotomous.csv data/phasmid_dyf data/MMPfiltered.SSID data/MMP_SNP_WeightFile.txt	


##======================================================================================
## Power Analysis
##======================================================================================

## Create Haplotype matrix and SNP location file




##======================================================================================
## Clean: files to delete if to reset to project start before analysis
##======================================================================================

clean:
	-rm -f data/phenotype_amphid_dyf_dichotomous.csv 
	-rm -f data/phenotype_phasmid_dyf_dichotomous.csv
	-rm -f data/list_VCstrains_vcf.txt
	-rm -f data/MMPcoding.vcf
	-rm -f data/MMPcoding_noHeader.vcf
	-rm -f data/MMPcoding.SSID
	## -rm -f data/MMPfiltered.vcf data/MMPfiltered.SSID
	-rm -f data/MMP_SNP_WeightFile.txt
	-rm -f data/amphid_dyf/MMPfiltered.fam data/amphid_dyf/MMPfiltered.bim data/amphid_dyf/MMPfiltered.bed data/amphid_dyf/MMPfiltered.log
	-rm -f data/amphid_dyf/SKAT_no_weights_results.txt data/amphid_dyf/SKAT_weights_results.txt data/amphid_dyf/SKAT_pANDq_no_weights_results.txt data/amphid_dyf/SKAT_pANDq_weights_results.txt
	-rmdir amphid/dyf
	-rm -f data/phasmid_dyf/MMPfiltered.fam data/phasmid_dyf/MMPfiltered.bim data/phasmid_dyf/MMPfiltered.bed data/phasmid_dyf/MMPfiltered.log
	-rm -f data/phasmid_dyf/SKAT_no_weights_results.txt data/phasmid_dyf/SKAT_weights_results.txt data/phasmid_dyf/SKAT_pANDq_no_weights_results.txt data/phasmid_dyf/SKAT_pANDq_weights_results.txt
	-rmdir phasmid/dyf