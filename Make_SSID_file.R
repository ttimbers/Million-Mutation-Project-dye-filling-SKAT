##Make snp(variant) set ID (SSID) file for SKAT

make_SSID <- function(vcf_sans_header) {
  ##takes merged .vcf file without header (arguement: vcf_sans_header) and creates a dataframe from which you 
  ##can save a snp(variant) set ID (SSID) file for SKAT analysis
  
  ##import vcf file without header
  mydata.vcf <- read.table(vcf_sans_header)

  ##get variant names
  variants <- mydata.vcf$V3

  ##get sequence names
  library(stringr)
  seq_tag  <- "SN=[a-zA-Z0-9.]{1,}"
  sequence_names  <- str_extract(mydata.vcf$V8, seq_tag)
  sequence_names  <- sub("SN=", "", sequence_names)

  vars <- as.character(variants)
  seqs  <- as.character(sequence_names)

  SSID  <- data.frame(cbind(seqs, vars))
  return(SSID)
}

##make SSID file for my set of MMP strains
##note - takes some time to read in MMPdyf_non-syn_coding_no_header.txt (derived from merged .vcf file)
my.SSID <- make_SSID("MMPdyf_non-syn_coding_no_header.txt")
##save SSID file
write.table(my.SSID, "MMPdyf_non-syn_coding_SSID.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, append=FALSE)

