# Tiffany Timbers, 2016
# Script to get proportion data
#
# Inputs:
#   1. a comma delimited .csv file (which has 5 columns: Strain,phenotype, N)
#   2. output filename
#
# Returns a tab delimited .tsv file (which has 4 columns: Strain, N, prop).

main <- function(){
  
  # assign command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  input_file_path <- args[1]
  output_phenotype_file <- args[2]
  
  # read in data
  phenotype_counts <- read.csv(input_file_path, header=TRUE)
  
  # remove VC2010
  phenotype_counts <- phenotype_counts[phenotype_counts$Strain != 'VC2010',]
  
  # add column of proportion of worms with defects
  phenotype_counts$prop <- phenotype_counts[,2]/phenotype_counts$N
  
  # remove unnecessary columns
  phenotype_counts[,2] <- NULL
  
  # save this data to a file to be used for SKAT analysis
  write.table(phenotype_counts, output_phenotype_file, sep="\t", row.names=FALSE, quote=FALSE, append=FALSE)
  
}

main()