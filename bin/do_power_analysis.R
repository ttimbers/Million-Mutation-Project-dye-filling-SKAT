## Tiffany Timbers
## August 3, 2015
##
## Power analysis for SKAT on the C. elegans Million Mutation Project Strains.
## Note - the median size of coding genes is 1,956 bp in C. elegans
## 
##
main <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  haplotype_file <- args[1]
  SNPlocation_file <- args[2]
  
  ## read in files
  haplotype <- read.table(haplotype_file, header=FALSE)
  SNPlocation <- read.table(SNPlocation_file, header=FALSE)
  
  ## Calculate the average power of randomly selected 1,956 bp regions
  ## with the following conditions.
  ##
  ## Subregion.Length = 1956
  ## Prevalence = 0.08
  ## Case.Prop = 0.08
  ## Causal percent = ?
  ## Causal.MAF.Cutoff = ?
  ## alpha = ?
  ## N.Sample.ALL = ?
  ## Weight.Param = ?
  ## N.Sim = 100
  ## OR.Type = ?
  ## MaxOR = ?
  ## Negative percent = ?
  ## r.corr = ?
 
  out.b <- Power_Logistic(Haplotypes = haplotype, SNP.Location = SNPlocation, SubRegion.Length = 1956, Prevalence = 0.08, Case.Prop=0.08)
  
}