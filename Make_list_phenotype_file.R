## Takes a file where the first column is strains, and subsequent columns are phenotypes and makes a file 
## called phenotypes.txt which lists the phenotypes (one per line)
##
## note - the input file containing the phenotypes must have a header which includes the names of the phenotypes
##
## Arguements: 
##    1. the name of the input file containing the phenotypes
##    2. the name of the output file
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  pheno.file <- args[1]
  pheno.out <- args[2]

  phenos <- read.table(pheno.file, header = TRUE)
  phenos.cols <- colnames(phenos)
  phenos.list <- phenos.cols[2:length(phenos.cols)]
  
  write(phenos.list,sep="\n", file=pheno.out)
}

main()