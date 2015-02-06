remove.vcf_header <- function(filename) {
##remove header from .vcf file
filename.out <- paste(substr(filename, 1, nchar(filename)-4), "_no_header.txt", sep="")
system(paste("grep -v '^#' ", filename, " > ", filename.out, sep=""))
}

remove.vcf_header("MMPdyf_non-syn_coding.vcf")
