# Million-Mutation-Project-dye-filling-SKAT
Genomic data and code to accompany the SKAT analysis of Million Mutation Project strains 
dye-filling phenotypes reported in Timbers et al.

All data used for this project is stored in a private repository on Figshare, which will 
be made public upon publication of the manuscript.

## Dependencies

`Bash Shell`, `Perl`, `Plink`, `R` and `R packages SKAT, stringr, fdrtool`

## How to use it

Call the Bash driver script, `bin/SKAT_MMP_dyf_driver.sh` and pass it the following 
arguments:

1. the path to the .vcf file containing all the genomic data for all the MMP strains 
	(e.g. MMP.vcf.gz )
2. the path to your phenotype file 
3. the path to the directory where you would like to save your results
4. Weight to be assigned to "Knockout" (e.g. nonsense, frameshift indels, splice site) 
	mutations
5. Weight to be assigned to in-frame indel mutations
6. Weight to be assigned to missense mutations

Example call:

`Bash bin/SKAT_MMP_dyf_driver.sh data data/MMP.vcf.gz data/dyf_phenotpe_with_N_and_wt.csv results/ 1 0.75 0.25