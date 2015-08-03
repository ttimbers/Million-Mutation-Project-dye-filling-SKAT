## Tiffany Timbers
## August 3, 2015
##
## Takes a vcf file, a SNP set ID file and the minimum number of variants you want to do SKAT for to create a vcf file and a a SNP set ID file with only those variants.

main <- function() {
  ## command line arguements
  args <- commandArgs(trailingOnly = TRUE)
  vcf_file <- args[1]
  SSID_file <- args[2]
  min_num_variants <- args[3]
  vcf_output_file <- args[4]
  SSID_output_file <- args[5]
  
  ## load files
  SSID <- read.table(SSID_file, header=FALSE)
  vcf <- read.table(vcf_file, header=FALSE)
  
  ## make a list of the variants to keep (i.e. count the number of occurrences of each gene)
  SSID_count <- plyr::count(SSID[,1])
  
  ## make SSID file with only those variants (i.e. reduce the SSID list to only
  ## those >= min_num_variants)
  SSID_count_greater_or_equal_to_min <- SSID[which(SSID_count[2] >= min_num_variants),]
  SSID_to_print <- data.frame(lapply(SSID_count_greater_or_equal_to_min, as.character), stringsAsFactors=FALSE)
  
  ## save SSID file
  write.table(SSID_to_print, file = SSID_output_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  ## make a .vcf file with only those variants
  require(dplyr)
  variant_list <- data.frame(SSID_count_greater_or_equal_to_min[,2])
  colnames(variant_list) <- ("V3")
  var_table <- tbl_df(variant_list)
  vcf_table <- tbl_df(vcf)
  vcf_count_greater_or_equal_to_min <- inner_join(vcf_table, var_table)
  vcf_to_print <- data.frame(lapply(vcf_count_greater_or_equal_to_min, as.character), stringsAsFactors=FALSE)
  
  ## save .vcf file
  system(paste("grep -E '^#' ", vcf_file, " > ", vcf_output_file, sep=""))
  write.table(vcf_count_greater_or_equal_to_min, file = vcf_output_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
}

main()