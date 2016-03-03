## Tiffany Timbers
## August 3, 2015
##
## Takes a vcf file, a SNP set ID file and the minimum number of variants you want to do SKAT for to create a vcf file and a a SNP set ID file with only those variants.

main <- function() {
  ## load libraries
  require(dplyr)

  ## command line arguements
  args <- commandArgs(trailingOnly = TRUE)
  vcf_file <- args[1]
  SSID_file <- args[2]
  # command line arg is a string, so need to be converted to a number
  min_num_variants <- as.numeric(args[3])
  vcf_output_file <- args[4]
  SSID_output_file <- args[5]

  ## load files
  SSID <- read.table(SSID_file, header=FALSE)
  colnames(SSID) <- c("gene", "variant")
  vcf <- read.table(vcf_file, header=FALSE)

  ## make a list of the variants to keep (i.e. count the number of occurrences of each gene)
  SSID_count <- plyr::count(SSID[,1])

  ## make SSID file with only those variants (i.e. reduce the SSID list to only
  ## those >= min_num_variants)
  genes_to_include <- data.frame(SSID_count[SSID_count$freq >= min_num_variants, 1])
  colnames(genes_to_include) <- ("gene")
  gene_table <- tbl_df(genes_to_include)
  SSID_table <- tbl_df(SSID)
  SSID_count_greater_or_equal_to_min <- semi_join(SSID_table, gene_table)
  SSID_to_print <- data.frame(lapply(SSID_count_greater_or_equal_to_min, as.character), stringsAsFactors=FALSE)

  ## save SSID file
  write.table(SSID_to_print, file = SSID_output_file, row.names = FALSE, col.names = FALSE, quote = FALSE)

  ## make a .vcf file with only those variants
  variant_list <- data.frame(SSID_count_greater_or_equal_to_min[,2])
  colnames(variant_list) <- ("V3")
  var_table <- tbl_df(variant_list)
  vcf_table <- tbl_df(vcf)
  vcf_count_greater_or_equal_to_min <- semi_join(vcf_table, var_table)
  vcf_to_print <- data.frame(lapply(vcf_count_greater_or_equal_to_min, as.character), stringsAsFactors=FALSE)

  ## save .vcf file
  system(paste("grep -E '^#' ", vcf_file, " > ", vcf_output_file, sep=""))
  write.table(vcf_count_greater_or_equal_to_min, file = vcf_output_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
}

main()

