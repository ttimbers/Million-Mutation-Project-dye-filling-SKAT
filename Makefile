## Tiffany A. Timbers
## July 25, 2015
##
## Makefile pipeline for SKAT analysis of genomic and phenotype data of Million Mutation
## Project strains dye-filling phenotypes reported in Timbers et al.
##
## Dependencies:
## 		Bash Shell (version 3.2.57(1))
## 		Make (version 3.81)
## 		Perl (v5.18.2)
## 		Plink (v1.90b3.36) (assumes that you have made plink executable and put it in your Bash Shell's $PATH)
## 		R (version 3.2.3)
## 		R packages:
## 		dplyr (version 0.4.3)
## 		fdrtool (version 1.2.15)
## 		plyr (version 1.8.3)
## 		pwr (version 1.1-3)
## 		SKAT (version 1.0.9)
## 		stringr (version 1.0.0)

## Final files to be generated from SKAT analysis
all: data/Table_S3.csv data/Table_S4.csv data/Table_S5.csv data/Table_S6.csv

##======================================================================================
## Analysis of amphid dye-filling phenotype
##======================================================================================

## Transforms data count data to a proportion. This script takes 2
## arguments: (1) a comma delimited .csv file (which has 3 columns: Strain,
## phenotype_count, N), and (2) output filename.
##
## Returns a tab delimited .tsv file (which has 3 columns: strain, N, prop).
data/phenotype_amphid_dyf_not_transformed.tsv: bin/get_proportions.R data/phenotype_amphid_dyf.csv
	Rscript bin/get_proportions.R data/phenotype_amphid_dyf.csv data/phenotype_amphid_dyf_not_transformed.tsv

## Transforms count data to a log transformed proportion. This script takes 2
## arguments: (1) a comma delimited .csv file (which has 3 columns: Strain,
## phenotype_count, N), and (2) output filename.
##
## Returns a tab delimited .tsv file (which has 4 columns: strain, N, log_prop, constant).
data/phenotype_amphid_dyf_log_05.tsv: bin/log_transform_phenotype.R data/phenotype_amphid_dyf.csv
	Rscript bin/log_transform_phenotype.R data/phenotype_amphid_dyf.csv data/phenotype_amphid_dyf_log_05.tsv 0.05 data/amphid_log_05

data/phenotype_amphid_dyf_log_005.tsv: bin/log_transform_phenotype.R data/phenotype_amphid_dyf.csv
	Rscript bin/log_transform_phenotype.R data/phenotype_amphid_dyf.csv data/phenotype_amphid_dyf_log_005.tsv 0.005 data/amphid_log_005


## Extracts the strain column from the data/phenotype_amphid_dyf_log.tsv file created
## above and saves it as list_VCstrains_vcf.txt
data/list_VCstrains_vcf.txt: data/phenotype_amphid_dyf_not_transformed.tsv
	awk '{print $$1}' data/phenotype_amphid_dyf_not_transformed.tsv | grep -h "^VC*" > data/list_VCstrains_vcf.txt

## Creates one merged .vcf file for only the strains that were assayed (listed in a file
## called list_VCstrains_vcf.txt. Uses -protein argument creates a merged .vcf file with
## only variants which affect coding regions
data/MMPcoding.vcf: data/MMP.vcf.gz bin/filter_MMP_variants.pl data/list_VCstrains_vcf.txt
	gunzip -c data/MMP.vcf.gz | perl bin/filter_MMP_variants.pl -input - -output data/MMPcoding.vcf -strain data/list_VCstrains_vcf.txt -protein

## Create SSID file for SKAT analysis
data/MMPcoding.SSID: bin/Make_SSID_file.R data/MMPcoding.vcf
	Rscript bin/Make_SSID_file.R data/MMPcoding.vcf data/MMPcoding.SSID

## Create a filtered SSID file and vcf file for only variants from those genes which have
## a specified minimum number of alleles (we chose 7)
data/MMPfiltered.vcf data/MMPfiltered.SSID: bin/create_reduced_variant_files.R data/MMPcoding.vcf data/MMPcoding.SSID
	Rscript bin/create_reduced_variant_files.R data/MMPcoding.vcf data/MMPcoding.SSID 7 data/MMPfiltered.vcf data/MMPfiltered.SSID

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
	if [ ! -d "data/amphid_dyf/" ]; then mkdir data/amphid_dyf; fi
	plink --vcf data/MMPfiltered.vcf --allow-extra-chr --no-fid --no-parents --no-sex --no-pheno --out data/amphid_dyf/MMPfiltered

## Perform linear regression SKAT analysis
data/amphid_dyf/SKAT_pANDq_no_weights_results.txt data/amphid_dyf/SKAT_pANDq_weights_results.txt: bin/do_linear_SKAT.R data/amphid_dyf/MMPfiltered.fam data/MMP_SNP_WeightFile.txt data/MMPfiltered.SSID data/phenotype_amphid_dyf_not_transformed.tsv
	Rscript bin/do_linear_SKAT.R data/amphid_dyf/MMPfiltered.fam data/phenotype_amphid_dyf_not_transformed.tsv data/amphid_dyf data/MMPfiltered.SSID data/MMP_SNP_WeightFile.txt

data/amphid_dyf/SKAT_pANDq_no_weights_results_log_05.txt data/amphid_dyf/SKAT_pANDq_weights_results_log_05.txt: bin/do_linear_SKAT.R data/amphid_dyf/MMPfiltered.fam data/MMP_SNP_WeightFile.txt data/MMPfiltered.SSID data/phenotype_amphid_dyf_log_05.tsv
	Rscript bin/do_linear_SKAT.R data/amphid_dyf/MMPfiltered.fam data/phenotype_amphid_dyf_log_05.tsv data/amphid_dyf data/MMPfiltered.SSID data/MMP_SNP_WeightFile.txt _log_05

data/amphid_dyf/SKAT_pANDq_no_weights_results_log_005.txt data/amphid_dyf/SKAT_pANDq_weights_results_log_005.txt: bin/do_linear_SKAT.R data/amphid_dyf/MMPfiltered.fam data/MMP_SNP_WeightFile.txt data/MMPfiltered.SSID data/phenotype_amphid_dyf_log_005.tsv
	Rscript bin/do_linear_SKAT.R data/amphid_dyf/MMPfiltered.fam data/phenotype_amphid_dyf_log_005.tsv data/amphid_dyf data/MMPfiltered.SSID data/MMP_SNP_WeightFile.txt _log_005

## Create Table S3 (Genome-wide association results from the SKAT of MMP DNA
## sequence variance and amphid dye-filling when variants were assigned biologically
## relevant weights. Results are sorted by p-value)
data/Table_S3.csv: bin/create_supp_results_table.R data/amphid_dyf/SKAT_pANDq_weights_results_log_005.txt data/gene_publicName_N_sequenceName.txt
	Rscript bin/create_supp_results_table.R data/amphid_dyf/SKAT_pANDq_weights_results_log_005.txt data/gene_publicName_N_sequenceName.txt data/Table_S3.csv

## Create Table S5 (Table S5. Genome-wide association results from the SKAT of MMP DNA
## sequence variance and amphid dye-filling when all variants were weighted equally.
## Results are sorted by p-value)
data/Table_S5.csv: bin/create_supp_results_table.R data/amphid_dyf/SKAT_pANDq_no_weights_results_log_005.txt data/gene_publicName_N_sequenceName.txt
	Rscript bin/create_supp_results_table.R data/amphid_dyf/SKAT_pANDq_no_weights_results_log_005.txt data/gene_publicName_N_sequenceName.txt data/Table_S5.csv


##======================================================================================
## Analysis of phasmid dye-filling phenotype
##======================================================================================

## Transforms data count data to a proportion. This script takes 2
## arguments: (1) a comma delimited .csv file (which has 3 columns: Strain,
## phenotype_count, N), and (2) output filename.
##
## Returns a tab delimited .tsv file (which has 3 columns: strain, N, prop).
data/phenotype_phasmid_dyf_not_transformed.tsv: bin/get_proportions.R data/phenotype_phasmid_dyf.csv
	Rscript bin/get_proportions.R data/phenotype_phasmid_dyf.csv data/phenotype_phasmid_dyf_not_transformed.tsv

## Assign dichotomous phasmid dye-filling phenotype to strains. This script takes 2
## arguments: (1) a comma delimited .csv file (which has 5 columns: Strain,
## phenotype, N), and (2) output filename.
##
## Returns a tab delimited .csv file (which has 4 columns: strain, phenotype, pvals,
## fdr.bh). The phenotype column contains 1 if strain's phenotype diverges significantly
## from wild-type (Fisher's exact test) and contain 1 if strain's phenotype diverges
## significantly from wild-type, and 0 if it does not. A 5% FDR (Benjamini-Hochberg
## procedure) is used to adjust for multiple comparisons.
data/phenotype_phasmid_dyf_log_05.tsv: bin/log_transform_phenotype.R data/phenotype_phasmid_dyf.csv
	Rscript bin/log_transform_phenotype.R data/phenotype_phasmid_dyf.csv data/phenotype_phasmid_dyf_log_05.tsv 0.05 data/phasmid_log_05

data/phenotype_phasmid_dyf_log_005.tsv: bin/log_transform_phenotype.R data/phenotype_phasmid_dyf.csv
	Rscript bin/log_transform_phenotype.R data/phenotype_phasmid_dyf.csv data/phenotype_phasmid_dyf_log_005.tsv 0.005 data/phasmid_log_005

## Create binary plink files for phasmid phenotype from filtered .vcf file
data/phasmid_dyf data/phasmid_dyf/MMPfiltered.fam data/phasmid_dyf/MMPfiltered.bim data/phasmid_dyf/MMPfiltered.bed data/phasmid_dyf/MMPfiltered.log: data/MMPfiltered.vcf
	if [ ! -d "data/phasmid_dyf/" ]; then mkdir data/phasmid_dyf; fi;
	plink --vcf data/MMPfiltered.vcf --allow-extra-chr --out data/phasmid_dyf/MMPfiltered

## Perform linear regression SKAT analysis
data/phasmid_dyf/SKAT_pANDq_no_weights_results.txt data/phasmid_dyf/SKAT_pANDq_weights_results.txt: bin/do_linear_SKAT.R data/phasmid_dyf/MMPfiltered.fam data/MMP_SNP_WeightFile.txt data/MMPfiltered.SSID data/phenotype_phasmid_dyf_not_transformed.tsv
	Rscript bin/do_linear_SKAT.R data/phasmid_dyf/MMPfiltered.fam data/phenotype_phasmid_dyf_not_transformed.tsv data/phasmid_dyf data/MMPfiltered.SSID data/MMP_SNP_WeightFile.txt

data/phasmid_dyf/SKAT_pANDq_no_weights_results_log_05.txt data/phasmid_dyf/SKAT_pANDq_weights_results_log_05.txt: bin/do_linear_SKAT.R data/phasmid_dyf/MMPfiltered.fam data/MMP_SNP_WeightFile.txt data/MMPfiltered.SSID data/phenotype_phasmid_dyf_log_05.tsv
	Rscript bin/do_linear_SKAT.R data/phasmid_dyf/MMPfiltered.fam data/phenotype_phasmid_dyf_log_05.tsv data/phasmid_dyf data/MMPfiltered.SSID data/MMP_SNP_WeightFile.txt _log_05

data/phasmid_dyf/SKAT_pANDq_no_weights_results_log_005.txt data/phasmid_dyf/SKAT_pANDq_weights_results_log_005.txt: bin/do_linear_SKAT.R data/phasmid_dyf/MMPfiltered.fam data/MMP_SNP_WeightFile.txt data/MMPfiltered.SSID data/phenotype_phasmid_dyf_log_005.tsv
	Rscript bin/do_linear_SKAT.R data/phasmid_dyf/MMPfiltered.fam data/phenotype_phasmid_dyf_log_005.tsv data/phasmid_dyf data/MMPfiltered.SSID data/MMP_SNP_WeightFile.txt _log_005


## Create Table S4 (Genome-wide association results from the SKAT of MMP DNA
## sequence variance and phasmid dye-filling when variants were assigned biologically
## relevant weights. Results are sorted by p-value)
data/Table_S4.csv: bin/create_supp_results_table.R data/phasmid_dyf/SKAT_pANDq_weights_results_log_05.txt data/gene_publicName_N_sequenceName.txt
	Rscript bin/create_supp_results_table.R data/phasmid_dyf/SKAT_pANDq_weights_results_log_05.txt data/gene_publicName_N_sequenceName.txt data/Table_S4.csv

## Create Table S6 (Genome-wide association results from the SKAT of MMP DNA
## sequence variance and phasmid dye-filling when all variants were weighted equally.
## Results are sorted by p-value)
data/Table_S6.csv: bin/create_supp_results_table.R data/phasmid_dyf/SKAT_pANDq_no_weights_results_log_05.txt data/gene_publicName_N_sequenceName.txt
	Rscript bin/create_supp_results_table.R data/phasmid_dyf/SKAT_pANDq_no_weights_results_log_05.txt data/gene_publicName_N_sequenceName.txt data/Table_S6.csv

##======================================================================================
## Bootstrap Power Analysis
##======================================================================================

## to run power analysis in a reasonable amount of time you need access to a cluster
## I used compute canada's Guillimin.

## The code I ran can be found in run_bootstrap_power.sh

##======================================================================================
## Clean: files to delete if to reset to project start before analysis
##======================================================================================

clean:
	# amphid associated files
	-rm -f data/phenotype_amphid_dyf_not_transformed.tsv
	-rm -f data/phenotype_amphid_dyf_log_05.tsv
	-rm -f data/phenotype_amphid_dyf_log_005.tsv
	-rm -f data/list_VCstrains_vcf.txt
	-rm -f data/MMPcoding.vcf
	-rm -f data/MMPcoding.SSID
	-rm -f data/MMPfiltered.vcf data/MMPfiltered.SSID
	-rm -f data/MMP_SNP_WeightFile.txt
	-rm -f data/amphid_dyf data/amphid_dyf/MMPfiltered.fam data/amphid_dyf/MMPfiltered.bim data/amphid_dyf/MMPfiltered.bed data/amphid_dyf/MMPfiltered.log
	-rm -f data/amphid_dyf/SKAT_no_weights_results.txt data/amphid_dyf/SKAT_weights_results.txt data/amphid_dyf/SKAT_pANDq_no_weights_results.txt data/amphid_dyf/SKAT_pANDq_weights_results.txt
	-rm -f data/amphid_dyf/file* data/MMPfiltered*
	-rm -f data/Table_S3.csv
	-rm -f data/Table_S5.csv
	-rmdir data/amphid_dyf

	# phasmid associated files
	-rm -f data/phenotype_phasmid_dyf_dichotomous.csv
	-rm -f data/phasmid_dyf/MMPfiltered.fam data/phasmid_dyf/MMPfiltered.bim data/phasmid_dyf/MMPfiltered.bed data/phasmid_dyf/MMPfiltered.log
	-rm -f data/phasmid_dyf/SKAT_no_weights_results.txt data/phasmid_dyf/SKAT_weights_results.txt data/phasmid_dyf/SKAT_pANDq_no_weights_results.txt data/phasmid_dyf/SKAT_pANDq_weights_results.txt
	-rm -f data/phasmid_dyf/file* Rplots.pdf
	-rm -f data/Table_S4.csv
	-rm -f data/Table_S6.csv
	-rmdir data/phasmid_dyf

	# bootstrap power analysis related files
	-rm -f data/temp_phenotype_amphid_dyf_log.tsv
	-rm -f data/temp_list_VCstrains_vcf.txt
	-rm -f data/temp_MMPfiltered.vcf
	-rm -f data/power data/power/temp_MMPfiltered.fam data/power/temp_MMPfiltered.bim data/power/temp_MMPfiltered.bed data/power/temp_MMPfiltered.log
	-rm -f data/power/temp_SKAT_pANDq_no_weights_results.csv
