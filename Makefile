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

## Perform SKAT analysis
data/amphid_dyf/SKAT_no_weights_results.txt data/amphid_dyf/SKAT_weights_results.txt data/amphid_dyf/SKAT_pANDq_no_weights_results.txt data/amphid_dyf/SKAT_pANDq_weights_results.txt: bin/do_SKAT.R data/amphid_dyf/MMPfiltered.fam data/MMP_SNP_WeightFile.txt data/MMPfiltered.SSID data/phenotype_amphid_dyf_dichotomous.csv
	Rscript bin/do_SKAT.R data/amphid_dyf/MMPfiltered.fam data/phenotype_amphid_dyf_dichotomous.csv data/amphid_dyf data/MMPfiltered.SSID data/MMP_SNP_WeightFile.txt

## Calculate effect size
data/amphid_dyf/MMPfiltered.gvsp data/amphid_dyf/MMPfiltered.effectsize: bin/calculate_effect_size.R data/phenotype_amphid_dyf_dichotomous.csv data/amphid_dyf/SKAT_pANDq_no_weights_results.txt data/MMPfiltered.SSID data/MMPfiltered.vcf
	Rscript bin/calculate_effect_size.R data/phenotype_amphid_dyf_dichotomous.csv data/amphid_dyf/SKAT_pANDq_no_weights_results.txt data/MMPfiltered.SSID data/MMPfiltered.vcf data/gene_publicName_N_sequenceName.txt data/amphid_dyf/MMPfiltered.gvsp data/amphid_dyf/MMPfiltered.effectsize

## Create Table S3 (Genome-wide association results from the SKAT of MMP DNA 
## sequence variance and amphid dye-filling when variants were assigned biologically 
## relevant weights. Results are sorted by p-value)
data/Table_S3.csv: bin/create_supp_results_table.R data/amphid_dyf/SKAT_weights_results.txt data/amphid_dyf/MMPfiltered.effectsize
	Rscript bin/create_supp_results_table.R data/amphid_dyf/SKAT_pANDq_weights_results.txt data/amphid_dyf/MMPfiltered.effectsize data/Table_S3.csv

## Create Table S5 (Table S5. Genome-wide association results from the SKAT of MMP DNA 
## sequence variance and amphid dye-filling when all variants were weighted equally. 
## Results are sorted by p-value)
data/Table_S5.csv: bin/create_supp_results_table.R data/amphid_dyf/SKAT_no_weights_results.txt data/amphid_dyf/MMPfiltered.effectsize
	Rscript bin/create_supp_results_table.R data/amphid_dyf/SKAT_pANDq_no_weights_results.txt data/amphid_dyf/MMPfiltered.effectsize data/Table_S5.csv


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
	if [ ! -d "data/phasmid_dyf/" ]; then mkdir data/phasmid_dyf; fi;
	plink --vcf data/MMPfiltered.vcf --allow-extra-chr --out data/phasmid_dyf/MMPfiltered

## Perform SKAT analysis (use SSID file made above)
data/phasmid_dyf/SKAT_no_weights_results.txt data/phasmid_dyf/SKAT_weights_results.txt data/phasmid_dyf/SKAT_pANDq_no_weights_results.txt data/phasmid_dyf/SKAT_pANDq_weights_results.txt: bin/do_SKAT.R data/phasmid_dyf/MMPfiltered.fam data/phenotype_phasmid_dyf_dichotomous.csv data/MMPfiltered.SSID data/MMP_SNP_WeightFile.txt
	Rscript bin/do_SKAT.R data/phasmid_dyf/MMPfiltered.fam data/phenotype_phasmid_dyf_dichotomous.csv data/phasmid_dyf data/MMPfiltered.SSID data/MMP_SNP_WeightFile.txt	

## Calculate effect size
data/phasmid_dyf/MMPfiltered.gvsp data/phasmid_dyf/MMPfiltered.effectsize: bin/calculate_effect_size.R data/phenotype_phasmid_dyf_dichotomous.csv data/phasmid_dyf/SKAT_pANDq_no_weights_results.txt data/MMPfiltered.SSID data/MMPfiltered.vcf
	Rscript bin/calculate_effect_size.R data/phenotype_phasmid_dyf_dichotomous.csv data/phasmid_dyf/SKAT_pANDq_no_weights_results.txt data/MMPfiltered.SSID data/MMPfiltered.vcf data/gene_publicName_N_sequenceName.txt data/phasmid_dyf/MMPfiltered.gvsp data/phasmid_dyf/MMPfiltered.effectsize

## Create Table S4 (Genome-wide association results from the SKAT of MMP DNA 
## sequence variance and phasmid dye-filling when variants were assigned biologically 
## relevant weights. Results are sorted by p-value)
data/Table_S4.csv: bin/create_supp_results_table.R data/phasmid_dyf/SKAT_weights_results.txt data/phasmid_dyf/MMPfiltered.effectsize
	Rscript bin/create_supp_results_table.R data/phasmid_dyf/SKAT_pANDq_weights_results.txt data/phasmid_dyf/MMPfiltered.effectsize data/Table_S4.csv

## Create Table S6 (Genome-wide association results from the SKAT of MMP DNA 
## sequence variance and phasmid dye-filling when all variants were weighted equally. 
## Results are sorted by p-value)
data/Table_S6.csv: bin/create_supp_results_table.R data/phasmid_dyf/SKAT_no_weights_results.txt data/phasmid_dyf/MMPfiltered.effectsize
	Rscript bin/create_supp_results_table.R data/phasmid_dyf/SKAT_pANDq_no_weights_results.txt data/phasmid_dyf/MMPfiltered.effectsize data/Table_S6.csv

##======================================================================================
## Bootstrap Power Analysis
##======================================================================================

## Do this for N = 50, 100, 200, 300, 400

## Make loop over these targets 1000 times

## Create list of randomly sampled strains (without replacement) & phenotype data from data/phenotype_amphid_dyf_dichotomous.csv
data/temp_phenotype_amphid_dyf_dichotomous.csv: bin/create_random_samples.R data/phenotype_amphid_dyf_dichotomous.csv
	Rscript bin/create_random_samples.R data/phenotype_amphid_dyf_dichotomous.csv \t TRUE 1 50 data/temp_phenotype_amphid_dyf_dichotomous.csv
	
## Create list of randomly selected strains from temp_phenotype_amphid_dyf.csv
#data/temp_list_VCstrains_vcf.txt: data/temp_phenotype_amphid_dyf_dichotomous.csv
#	awk '{print $$1}' data/temp_phenotype_amphid_dyf_dichotomous.csv | grep -h "^VC*" > data/temp_list_VCstrains_vcf.txt

## Create a vcf file from these random selected strains, only variants from data/MMPfiltered.vcf
#data/temp_MMPfiltered.vcf: bin/filter_MMP_variants.pl data/MMPfiltered.vcf
#	gunzip -c data/MMP.vcf.gz | perl bin/filter_MMP_variants.pl -input - -output data/temp_MMPcoding.vcf -strain data/temp_list_VCstrains_vcf.txt -protein

## Create binary plink files from the vcf file
#data/power data/power/temp_MMPfiltered.fam data/power/temp_MMPfiltered.bim data/power/temp_MMPfiltered.bed data/power/temp_MMPfiltered.log: data/temp_MMPfiltered.vcf
#	if [ ! -d "data/power/" ]; then mkdir data/power; fi
#	plink --vcf data/temp_MMPfiltered.vcf --allow-extra-chr --no-fid --no-parents --no-sex --no-pheno --out data/power/temp_MMPfiltered

## Perform SKAT analysis
#data/power/SKAT_pANDq_no_weights_results.csv: bin/do_SKAT_no_weights.R data/power/temp_MMPfiltered.fam data/MMPfiltered.SSID data/phenotype_amphid_dyf_dichotomous.csv
#	Rscript bin/do_SKAT_no_weights.R data/amphid_dyf/MMPfiltered.fam data/phenotype_amphid_dyf_dichotomous.csv data/amphid_dyf data/MMPfiltered.SSID data/MMP_SNP_WeightFile.txt

##======================================================================================
## Plot data for paper (characterization of bgnt-1)
##======================================================================================
## Note - these are not in make syntax yet...

## plot ADL cilia length in bgnt-1 and wild-type
#Rscript bin/analyze_length.R data/ADL_cilia_length_unblind_2015-07-31.csv cilia data/cilia_length.pdf data/cilia_length.stats

## plot distal tip of ADL cilia to distal end of socket cell length in bgnt-1 and wild-type
#Rscript bin/analyze_length.R data/distal_end_of_ADL_cilia_to_distal_end_of_socket_unblind.csv cilia_socket_diss data/distal_end_of_ADL_cilia_to_distal_end_of_socket.pdf data/distal_end_of_ADL_cilia_to_distal_end_of_socket.stats

## plot proportion of animals that have > 1 ADL cilia per amphid
#Rscript bin/ADL-guidance.R data/ADL-guidance.csv data/ADL-guidance.pdf data/ADL-guidance.stats

##======================================================================================
## Clean: files to delete if to reset to project start before analysis
##======================================================================================

clean:
	# amphid associated files
	-rm -f data/phenotype_amphid_dyf_dichotomous.csv 
	-rm -f data/list_VCstrains_vcf.txt
	-rm -f data/MMPcoding.vcf
	-rm -f data/MMPcoding.SSID
	-rm -f data/MMPfiltered.vcf data/MMPfiltered.SSID
	-rm -f data/MMP_SNP_WeightFile.txt
	-rm -f data/amphid_dyf data/amphid_dyf/MMPfiltered.fam data/amphid_dyf/MMPfiltered.bim data/amphid_dyf/MMPfiltered.bed data/amphid_dyf/MMPfiltered.log
	-rm -f data/amphid_dyf/SKAT_no_weights_results.txt data/amphid_dyf/SKAT_weights_results.txt data/amphid_dyf/SKAT_pANDq_no_weights_results.txt data/amphid_dyf/SKAT_pANDq_weights_results.txt
	-rm -f data/amphid_dyf/file* data/MMPfiltered*
	-rm -f data/amphid_dyf/MMPfiltered.gvsp data/amphid_dyf/MMPfiltered.effectsize
	-rm -f data/Table_S3.csv
	-rm -f data/Table_S5.csv
	-rmdir data/amphid_dyf
	
	# phasmid associated files
	-rm -f data/phenotype_phasmid_dyf_dichotomous.csv
	-rm -f data/phasmid_dyf/MMPfiltered.fam data/phasmid_dyf/MMPfiltered.bim data/phasmid_dyf/MMPfiltered.bed data/phasmid_dyf/MMPfiltered.log
	-rm -f data/phasmid_dyf/SKAT_no_weights_results.txt data/phasmid_dyf/SKAT_weights_results.txt data/phasmid_dyf/SKAT_pANDq_no_weights_results.txt data/phasmid_dyf/SKAT_pANDq_weights_results.txt
	-rm -f data/phasmid_dyf/file* Rplots.pdf
	-rm -f data/phasmid_dyf/MMPfiltered.gvsp data/phasmid_dyf/MMPfiltered.effectsize
	-rm -f data/Table_S4.csv
	-rm -f data/Table_S6.csv
	-rmdir data/phasmid_dyf 
	
	# bootstrap power analysis related files
	-rm -f data/temp_phenotype_amphid_dyf_dichotomous.csv