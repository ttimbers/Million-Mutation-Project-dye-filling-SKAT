## Tiffany Timbers
## August 3, 2015
##
## Power analysis for SKAT on the C. elegans Million Mutation Project Strains.
## Note - the median size of coding genes is 1,956 bp in C. elegans
## 
##
## use below when run inside RStudio:
## haplotype_file <- "data/haplotype.matrix"
## SNPlocation_file <- "data/SNPlocation.file"
main <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  haplotype_file <- args[1]
  SNPlocation_file <- args[2]
  
  library(SKAT)
  
  ## read in files
  haplotype <- as.matrix(read.table(haplotype_file, header=FALSE))
  SNPlocation <- read.table(SNPlocation_file, header=FALSE)
  
  #data("SKAT.haplotypes")
  #haplotype_test <- SKAT.haplotypes$Haplotype[1:1000,]
  #SNPlocation_test <- SKAT.haplotypes$SNPInfo[,3]
  
  #haplotype <- haplotype[1:2000,]
  
  
  ## Calculate the average power of randomly selected 1,956 bp regions
  ## with the following conditions.
  ##
  ## Subregion.Length = 1956
  ## Prevalence = 0.08
  ## Case.Prop = 0.08
  ## Causal percent = 40
  ## Causal.MAF.Cutoff = 0.05
  ## alpha = 0.002
  ## N.Sample.ALL = 480
  ## Weight.Param = c(1,25)
  ## N.Sim = 100
  ## OR.Type = "Log"
  ## MaxOR = 13
  ## Negative percent = 0
  
  out.b <- Power_Logistic(Haplotypes = haplotype, SNP.Location = SNPlocation, SubRegion.Length = 1956, Prevalence = 0.08, Case.Prop=0.08, Causal.Percent = 40, Causal.MAF.Cutoff = 0.05, alpha = 0.002, Weight.Param = c(1,25), N.Sim = 100, OR.Type = "Log", MaxOR = 13, Negative.Percent = 0)
  
  out.b<-Power_Logistic(Haplotypes = haplotype, SNP.Location = SNPlocation, SubRegion.Length=1956, Case.Prop=0.08, Causal.Percent= 40, N.Sim=100, Causal.MAF.Cutoff = 0.05, MaxOR=13,Negative.Percent=0)
  
  Get_RequiredSampleSize(out.b, Power=0.5)
}

main()