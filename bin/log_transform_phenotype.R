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
  
  #load libraries
  library(ggplot2)
  
  # assign command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  input_file_path <- args[1]
  output_phenotype_file <- args[2]
  output_histogram_prefix <- args[3]
  
  # read in data
  phenotype_counts <- read.csv(input_file_path, header=TRUE)

  # remove VC2010
  phenotype_counts <- phenotype_counts[phenotype_counts$Strain != 'VC2010',]
  
  # add column of proportion of worms with defects
  phenotype_counts$prop <- phenotype_counts[,2]/phenotype_counts$N
  
  # make histogram of raw proportion data and save it
  raw_hist <- ggplot(data=phenotype_counts, aes(phenotype_counts$prop)) + 
    geom_histogram(breaks=seq(0, 1, by = 0.05), colour="black", fill="white") +
    labs(x="Proportion of dye-filling defects", y="Count")
  ggsave(filename = paste0(output_histogram_prefix, '_raw_phenotype_histogram.pdf'), height = 3.5, width = 3.5)
  
  # add 0.05 to each proportion so that there are no 0's
  phenotype_counts$log_prop <- phenotype_counts$prop + 0.05
  
  # add column of log transformed proportions
  phenotype_counts$log_prop <- log(phenotype_counts$log_prop)
  
  #phenotype_counts$log_prop <- asin(sqrt(phenotype_counts$log_prop))
  
  # make histogram of raw proportion data and save it
  log_hist <- ggplot(data=phenotype_counts, aes(phenotype_counts$log_prop)) + 
    geom_histogram(breaks=seq(-2.5, 0.5, by = 0.15), colour="black", fill="white") +
    labs(x="log(Proportion of dye-filling defects + 5e-03)", y="Count")
  log_hist
  ggsave(filename = paste0(output_histogram_prefix, '_log_transformed_phenotype_histogram.pdf'), height = 3.5, width = 3.5)

  # remove unnecessary columns
  phenotype_counts[,2] <- NULL
  phenotype_counts$N <- NULL
  
  # save this data to a file to be used for SKAT analysis
  write.table(phenotype_counts, output_phenotype_file, sep="\t", row.names=FALSE, quote=FALSE, append=FALSE)
        
}

main()