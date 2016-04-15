## Tiffany Timbers
## August 3, 2015
##
## Power analysis for SKAT on the C. elegans Million Mutation Project Strains.
## 
## 
## Approach - bootstap
## 0. Make a compressed vcf with only strains & snps that I have used in my analysis
##
## 1. randomly sample 50 worms strains from column 1 of data/phenotype_amphid_dyf.csv
##		* R script: randomly sample (without replacement) from a list data/phenotype_amphid_dyf.csv, save as temp file (with unique ID?)
##		* make temp vcf with these samples
## 		* convert files with plink and add phenotypes

## 2. do SKAT on amphid data, get p & q values, write those with q's < 0.3 to file naming it with time-date_stamp + sample N (here 100)
## 3. do this 10000 times
## 4. Report What proportion of the time do you get a significant q-value under 30%, under 20%, under 10%?
## 5. Do this all again for 100, 200 & 300 samples

