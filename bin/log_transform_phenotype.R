# Tiffany Timbers, 2016
# Script to log transform proportion data
#
# Inputs:
#   1. a comma delimited .csv file (which has 5 columns: Strain,phenotype, N)
#   2. output filename
#
# Returns a tab delimited .tsv file (which has 4 columns: strain, phenotype, dyf_proportion, 
# and logit(dyf_proportion)).

main <- function(){
  
  # assign command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  input_file_path <- args[1]
  output_file_path <- args[2]
  
  # read in data
  phenotype_counts <- read.csv(input_file_path, header=TRUE)

  # remove VC2010
  phenotype_counts <- phenotype_counts[phenotype_counts$Strain != 'VC2010',]
  
  # add column of proportion of worms with defects
  phenotype_counts$prop <- phenotype_counts$Amphid_defect_count/phenotype_counts$N
  
  # add 0.05 to each proportion so that there are no 0's
  phenotype_counts$prop <- phenotype_counts$prop + 0.05

  # add column of logit transformed proportions
  phenotype_counts$log_prop <- log(phenotype_counts$prop)

  # remove unnecessary columns
  phenotype_counts$Amphid_defect_count <- NULL
  phenotype_counts$N <- NULL
  
  # save this data to a file to be used for SKAT analysis
  write.table(phenotype_counts, output_file_path, sep="\t", row.names=FALSE, quote=FALSE, append=FALSE)
        
}

main()