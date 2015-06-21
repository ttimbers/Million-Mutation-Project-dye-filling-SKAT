all: data/phenotype_amphid_dyf_dichotomous.csv

data/phenotype_amphid_dyf_dichotomous.csv: bin/Assign_dichotomous_phenotype.R
	Rscript bin/Assign_dichotomous_phenotype.R data/phenotype_amphid_dyf.csv data/phenotype_amphid_dyf_dichotomous.csv
	
clean:
	-rm -f data/phenotype_amphid_dyf_dichotomous.csv
