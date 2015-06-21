all: data/filteredMMP.vcf

data/phenotype_amphid_dyf_dichotomous.csv: bin/Assign_dichotomous_phenotype.R
	Rscript bin/Assign_dichotomous_phenotype.R data/phenotype_amphid_dyf.csv data/phenotype_amphid_dyf_dichotomous.csv

data/list_VCstrains_vcf.txt: data/phenotype_amphid_dyf_dichotomous.csv
	awk '{print $$1}' data/phenotype_amphid_dyf_dichotomous.csv | grep -h "^VC*" > data/list_VCstrains_vcf.txt
	
data/filteredMMP.vcf: data/MMP.vcf.gz bin/filter_MMP_variants.pl data/list_VCstrains_vcf.txt	
	gzcat data/MMP.vcf.gz | perl bin/filter_MMP_variants.pl -input - -output data/filteredMMP.vcf -strain data/list_VCstrains_vcf.txt -protein
	
clean:
	-rm -f data/phenotype_amphid_dyf_dichotomous.csv data/filteredMMP.vcf