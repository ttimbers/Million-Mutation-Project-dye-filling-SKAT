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

data/phenotype_amphid_dyf_dichotomous.csv: bin/Assign_dichotomous_phenotype.R data/phenotype_amphid_dyf.csv
	Rscript bin/Assign_dichotomous_phenotype.R data/phenotype_amphid_dyf.csv data/phenotype_amphid_dyf_dichotomous.csv

data/list_VCstrains_vcf.txt: data/phenotype_amphid_dyf_dichotomous.csv
	awk '{print $$1}' data/phenotype_amphid_dyf_dichotomous.csv | grep -h "^VC*" > data/list_VCstrains_vcf.txt
	
data/filteredMMP.vcf: data/MMP.vcf.gz bin/filter_MMP_variants.pl data/list_VCstrains_vcf.txt	
	gzcat data/MMP.vcf.gz | perl bin/filter_MMP_variants.pl -input - -output data/filteredMMP.vcf -strain data/list_VCstrains_vcf.txt -protein
	
data/filteredMMP_noHeader.vcf: data/filteredMMP.vcf
	grep -v '^#' data/filteredMMP.vcf > data/filteredMMP_noHeader.vcf	
	
data/MMP_SNP_WeightFile.txt: bin/Make_custom_variant_weights.R data/filteredMMP_noHeader.vcf
	Rscript bin/Make_custom_variant_weights.R data/filteredMMP_noHeader.vcf 1 0.75 0.25 data/MMP_SNP_WeightFile.txt

data/amphid_dyf/filteredMMP.fam data/amphid_dyf/filteredMMP.bim data/amphid_dyf/filteredMMP.bed data/amphid_dyf/filteredMMP.log: data/filteredMMP.vcf
	mkdir data/amphid_dyf
	plink --vcf data/filteredMMP.vcf --allow-extra-chr --out data/amphid_dyf/filteredMMP

data/MMP_non-syn_coding_SSID.txt: bin/Make_SSID_file.R data/filteredMMP_noHeader.vcf
	Rscript bin/Make_SSID_file.R data/filteredMMP_noHeader.vcf data/MMP_non-syn_coding_SSID.txt

data/amphid_dyf/SKAT_no_weights_results.txt data/amphid_dyf/SKAT_weights_results.txt data/amphid_dyf/SKAT_pANDq_no_weights_results.txt data/amphid_dyf/SKAT_pANDq_weights_results.txt: bin/do_SKAT.R data/amphid_dyf/filteredMMP.fam data/MMP_SNP_WeightFile.txt data/MMP_non-syn_coding_SSID.txt
	Rscript bin/do_SKAT.R data/amphid_dyf/filteredMMP.fam data/phenotype_amphid_dyf_dichotomous.csv data/amphid_dyf data/MMP_non-syn_coding_SSID.txt data/MMP_SNP_WeightFile.txt

##======================================================================================
## Analysis of phasmid dye-filling phenotype
##======================================================================================

data/phenotype_phasmid_dyf_dichotomous.csv: bin/Assign_dichotomous_phenotype.R data/phenotype_phasmid_dyf.csv
	Rscript bin/Assign_dichotomous_phenotype.R data/phenotype_phasmid_dyf.csv data/phenotype_phasmid_dyf_dichotomous.csv

data/phasmid_dyf/filteredMMP.fam data/phasmid_dyf/filteredMMP.bim data/phasmid_dyf/filteredMMP.bed data/phasmid_dyf/filteredMMP.log: data/filteredMMP.vcf
	mkdir data/phasmid_dyf
	plink --vcf data/filteredMMP.vcf --allow-extra-chr --out data/phasmid_dyf/filteredMMP

data/phasmid_dyf/SKAT_no_weights_results.txt data/phasmid_dyf/SKAT_weights_results.txt data/phasmid_dyf/SKAT_pANDq_no_weights_results.txt data/phasmid_dyf/SKAT_pANDq_weights_results.txt: bin/do_SKAT.R data/phasmid_dyf/filteredMMP.fam data/MMP_SNP_WeightFile.txt data/MMP_non-syn_coding_SSID.txt
	Rscript bin/do_SKAT.R data/phasmid_dyf/filteredMMP.fam data/phenotype_phasmid_dyf_dichotomous.csv data/phasmid_dyf data/MMP_non-syn_coding_SSID.txt data/MMP_SNP_WeightFile.txt	


## files to delete if to reset to project start before analysis
clean:
	-rm -f data/phenotype_amphid_dyf_dichotomous.csv data/list_VCstrains_vcf.txt data/filteredMMP.vcf
	-rm -f data/filteredMMP_noHeader.vcf data/MMP_SNP_WeightFile.txt
	-rm -f data/amphid_dyf/filteredMMP.fam data/amphid_dyf/filteredMMP.bim data/amphid_dyf/filteredMMP.bed data/amphid_dyf/filteredMMP.log
	-rm -f data/amphid_dyf/file*
	-rmdir data/amphid_dyf
	-rm -f data/MMP_non-syn_coding_SSID.txt
	-rm -f data/amphid_dyf/SKAT_no_weights_results.txt data/amphid_dyf/SKAT_weights_results.txt data/amphid_dyf/SKAT_pANDq_no_weights_results.txt data/amphid_dyf/SKAT_pANDq_weights_results.txt
	-rm -f data/MMP_non-syn_coding_SSID.txt.temp 
	-rm -f data/filteredMMP_noHeader.coding data/filteredMMP_noHeader.indels data/filteredMMP_noHeader.snvs
	-rm -f data/phasmid_dyf/filteredMMP.fam data/phasmid_dyf/filteredMMP.bim data/phasmid_dyf/filteredMMP.bed data/phasmid_dyf/filteredMMP.log
	-rm -f data/phasmid_dyf/file*
	-rmdir data/phasmid_dyf
	-rm -f data/phasmid_dyf/SKAT_no_weights_results.txt data/phasmid_dyf/SKAT_weights_results.txt data/phasmid_dyf/SKAT_pANDq_no_weights_results.txt data/phasmid_dyf/SKAT_pANDq_weights_results.txt