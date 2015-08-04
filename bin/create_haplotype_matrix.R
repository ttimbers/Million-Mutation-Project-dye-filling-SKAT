## Tiffany Timbers
## July 29, 2015
##
## Takes a variant call file and transforms it into a haplo-type matrix to be used for SKAT power
## analysis. A haplotype matrix with each row as a different individual and each column as a 
## separate SNP (default= NULL). Each element of the matrix should be either 0 (major allel) or 
## 1 (minor allele). 
##
## Arguements:  (1) path to vcf file
##              (2) column number of the genotype for the first individual
##              (3) path and name of file to save Haplotype matrix to
##              (4) path and name of file to save SNP location vector to
##
## Dependencies:  (1) R and R libraries stringr and plyr
##                (2) Bash Shell

main <- function() {
  ## command line arguements
  args <- commandArgs(trailingOnly = TRUE)
  vcf_file <- args[1]
  haplotype_output_file <- args[2]
  SNP_location_output_file <- args[3]
  
  ## load libraries
  require(stringr)
  require(plyr)
  
  ## Make haplotype matrix
  
  ## read rectangular version into R
  vcf <- read.table(vcf_file, header=FALSE)
  
  ## reduce matrix
  temp_df <- df <- subset(vcf, select = -c(V1, V2, V4, V5, V6, V7, V8, V9))
  
  ## replace 0/0 with 0, and 1/1 with 1
  temp_df <- data.frame(lapply(temp_df, as.character), stringsAsFactors=FALSE)
  temp_df[temp_df=="0/0"] <- 0
  temp_df[temp_df=="1/1"] <- 1
  
  ## make a dataframe of the 0 and 1 values (i.e. ID column)
  haploytpe_matrix_t <- temp_df[,-1]
  
  ## transpose dataframe
  haplotype_matrix <- as.data.frame(t(as.matrix(haploytpe_matrix_t)))
  
  ## save to haplotype matrix to output file
  write.table(haplotype_matrix, file = haplotype_output_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  ## Make SNP location vector
  SNP_location_vector <- vcf[,2]
  
  ## save to SNP location vector to output file
  write.table(SNP_location_vector, file = SNP_location_output_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
}

main()