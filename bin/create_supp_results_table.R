## Tiffany Timbers
## August 7, 2015
##
## This script takes in the SKAT results from do_SKAT.R and the effect size calculated with 
## calculate_effect_size.R and creates the supplemental tables for the paper containing columns
## gene(seq name), p-value, q-value, number of variants, number of dyf strains with mutations in 
## that gene, number of non-dyf strains with mutations in that gene, and the effect size of that gene.

main <- function(){
  ## Get command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  SKAT_results_file <- args[1] 
  effectsize_file <- args[2] 
  outout_table_name <- args[3] 
  
  ## load libraries
  require(dplyr)
  
  ## load data files
  SKAT_results <- read.table(SKAT_results_file, header = TRUE)
  effectsize <- read.table(effectsize_file, header = TRUE)
  
  ## drop unnecessary columns
  SKAT_results$N.Marker.Test <- NULL
  SKAT_results$N.Marker.All <- NULL
  effectsize_tmp <- effectsize[,1:5]
  effectsize_tmp$effect_size <- effectsize[,10]
  effectsize <- effectsize_tmp
  rm(effectsize_tmp)
  
  ## rename columns
  colnames(SKAT_results) <- c("sequence", "p-value", "q-value")
  colnames(effectsize) <- c("public_gene_name", "sequence", "number_of_strains_with_mutations", "number_of_strains_with_dyf_phenotype", "number_of_strains_with_wild-type_phenotype", "effects_size")
  
  ## join SKAT_results and effectsize
  joined_data <- left_join(SKAT_results, effectsize, by = "sequence")
  
  ## order columns
  joined_data <- joined_data[c("public_gene_name", "sequence", "p-value", "q-value", "effects_size", "number_of_strains_with_mutations", "number_of_strains_with_dyf_phenotype", "number_of_strains_with_wild-type_phenotype")]
  
  ## write out table
  write.table(joined_data, outout_table_name, row.names = FALSE, quote = FALSE, append = FALSE)
}

main()