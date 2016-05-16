## Tiffany Timbers
## August 7, 2015
##
## This script takes in the SKAT results from do_SKAT.R and the effect size calculated with 
## calculate_effect_size.R and creates the supplemental tables for the paper containing columns
## gene(seq name), p-value, Bonferroni-adjusted-p-value, number of variants, number of dyf strains 
## with mutations in that gene, number of non-dyf strains with mutations in that gene

main <- function(){

  ## Start the clock!
  ptm <- proc.time()

  ## load libraries
  require(dplyr)

  ## Get command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  SKAT_results_file <- args[1] 
  gene_publicName_N_sequenceName <- args[2]
  outout_table_name <- args[3] 
  
  ## load data files
  SKAT_results <- read.table(SKAT_results_file, header = TRUE)
  gene_publicName_N_sequenceName <- read.table(gene_publicName_N_sequenceName, header=FALSE, sep="\t")
  
  ## drop unnecessary columns
  SKAT_results$N.Marker.All <- NULL
  
  ## rename columns
  colnames(SKAT_results) <- c("sequence", "p-value", "number_of_strains_with_mutations", "Bonferroni_adjusted_p-value")
  colnames(gene_publicName_N_sequenceName) <- c("public_gene_name", "sequence")
  
  ## join tables
  joined_data <- left_join(SKAT_results, gene_publicName_N_sequenceName, by = "sequence")
  
  ## order columns
  joined_data <- joined_data[c("public_gene_name", "sequence", "p-value", "Bonferroni_adjusted_p-value", "number_of_strains_with_mutations")]
  
  ## write out table
  write.table(joined_data, outout_table_name, row.names = FALSE, quote = FALSE, append = FALSE)
  
  ## output time to run script
  the_time <- proc.time() - ptm # Stop the clock
  print(paste("It took", the_time[3], "to run create_supp_table.R"))
}

main()