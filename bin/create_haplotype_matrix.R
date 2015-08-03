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
##              (3) path and name of file to save to
##
## Dependencies:  (1) R and R libraries stringr and plyr
##                (2) Bash Shell

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  vcf_file <- args[1]
  first_genotype_column <- args[2]
  output_file <- args[3]
  
  ## load libraries
  require(stringr)
  require(plyr)
  
  ## reduce matrix to only variant names and variant calls
  
  ## make a variable to hold the filename I will give to a rectangular version of the vcf 
  ## file (withoug column 1)
  vcf_df_no_header_name <- paste(vcf_file, ".no_header", sep="")
  
  ## make a rectangular version of the vcf file (withoug column 1)
  system(paste("grep -E -v '^##' ", vcf_file, " | cut -d ' ' -f 2,3,", first_genotype_column, 
               "- > ", vcf_df_no_header_name, sep = ''))
  
  ## read rectangular version into R
  vcf_df_no_header <- read.table(vcf_df_no_header_name, header=TRUE)
  
  ## reduce matrix
  temp_df <- vcf_df_no_header[2:dim(vcf_df_no_header)[2]]
  
  ## replace 0/0 with 0, and 1/1 with 1
  temp_df <- data.frame(lapply(temp_df, as.character), stringsAsFactors=FALSE)
  temp_df[temp_df=="0/0"] <- 0
  temp_df[temp_df=="1/1"] <- 1
  
  ## make column #1 the row names
  temp_df2 <- temp_df[,-1]
  rownames(temp_df2) <- temp_df[,1]
  
  # transpose dataframe
  haplotype_matrix <- as.data.frame(t(as.matrix(temp_df2[-1])))
  
  # save to output file
  write.table(haplotype_matrix, file = output_file, row.names = TRUE, col.names = TRUE, 
              quote = FALSE)
}

main()