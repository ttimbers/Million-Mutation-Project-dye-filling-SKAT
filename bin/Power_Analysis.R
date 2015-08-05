## Tiffany Timbers, July 5, 2015
## Power analysis of SKAT for dye-filling phenotypes from 485 MMP strains from 
## Timbers et al. 2005 

## Depends on a .vcf file with only a oneline header that contains strain names
## This file was made with the following command in the Shell:
## grep -v "##" data/filteredMMP.vcf > data/filteredMMP_oneline_header.vcf
## Also needed to delete first character of file, did that using the following unix 
## commands:
## head -1 data/filteredMMP_oneline_header.vcf | awk '{ gsub("#", "", $1); print $0}' > data/MMP_oneline_header.vcf
## grep -v "#" data/filteredMMP_oneline_header.vcf >> data/MMP_oneline_header.vcf
## rm data/filteredMMP_oneline_header.vcf

## use below when run inside RStudio:
setwd("Documents/Post-Doc/Manuscripts/MMP_dyf_screen/code/Million-Mutation-Project-dye-filling-SKAT/")
path_to_phenotypes <- "data/phenotype_amphid_dyf_dichotomous.csv"
path_to_SKAT_analysis <- "data/amphid_dyf/SKAT_version1_0_9_results/SKAT_pANDq_no_weights_results.txt"
path_to_SSID <- "data/MMPfiltered.SSID"
path_to_vcf <- "data/MMPfiltered.vcf"
output_gvsp_file <- "data/MMPfiltered.gvsp"
output_effect_size_file <- "data/MMPfiltered.effectsize"

main <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  path_to_phenotypes <- args[1]
  path_to_SKAT_analysis <- args[2]
  path_to_SSID <- args[3]
  path_to_vcf <- args[4]
  output_gvsp_file <- args[5]
  output_effect_size_file <- args[6]
  
  require(pwr)
  require(dplyr)
  require(plyr)
  
  ## Read in phenotype file
  phenotypes <- read.csv(path_to_phenotypes, header = TRUE, sep = "\t")
  
  ## Read in SKAT analysis results file
  SKAT_results <- read.table(path_to_SKAT_analysis, header = TRUE, sep="\t")
  
  ## Read in SSID file
  SSID <- read.table(path_to_SSID)
  colnames(SSID) <- c("gene", "variant")
  
  ## Read in VCF file with oneline header
  ## make a variable to hold the filename I will give to a rectangular version of the vcf 
  ## file (withoug column 1)
  vcf_df_no_header_name <- paste(path_to_vcf, ".no_header", sep="")
 
  ## make a rectangular version of the vcf file (withoug column 1)
  system(paste("grep -E -v '^##' ", path_to_vcf, " | cut -d$'\t' -f 3,10- > ", vcf_df_no_header_name, sep = ''))
  
  vcf <- read.table(vcf_df_no_header_name, header=TRUE)
  
  ## make variants row names
  rownames(vcf) <- vcf[,1]
  vcf[,1] <- NULL
  
  ## Change from factors to character (and preserve row names)
  i <- sapply(vcf, is.factor)
  vcf[i] <- lapply(vcf[i], as.character)
  
  ## get row and column names of all 1/1 values
  var_and_strain <- Which.names(vcf, value="1/1")
  colnames(var_and_strain) <- c("variant", "strain")
  
  ## Add genes to data frame
  gvsp <- left_join(var_and_strain, SSID)
  gvsp <- gvsp[c("gene", "variant", "strain")]
  
  ## Add phenotypes to data frame
  ps <- phenotypes[,1:2]
  gvsp <- left_join(gvsp, phenotypes[,1:2])
  gvsp <- gvsp[c("gene", "variant", "strain")]
  
  ## save gvsp file
  write.table(gvsp, output_gvsp_file, row.names = FALSE, quote = FALSE, append = FALSE)
  
  ## Determine average effect size for all genes variants from SKAT analysis
  effect_size <- ddply(gvsp, c("gene"), summarise,
                 N_variants  = length(variant),
                 N_cases = sum(phenotype))
  
  effect_size$N_controls <- effect_size$N_variants - effect_size$N_cases
  effect_size$prob_minor <- effect_size$N_cases/effect_size$N_variants
  effect_size$prob_major <- (40 - effect_size$N_cases) / ((40 - effect_size$N_cases) + (440 - effect_size$N_controls))
  effect_size$odds_minor <- effect_size$prob_minor / (1 - effect_size$prob_minor)
  effect_size$odds_major <- effect_size$prob_major / (1 - effect_size$prob_major)
  effect_size$odds_ratio <- effect_size$odds_minor / effect_size$odds_major 
  
  write.table(effect_size, output_effect_size_file, row.names = FALSE, quote = FALSE, append = FALSE)
  
  ## Make a list of genes (& their number of vars) that are significantly associated 
  ## with the phenotype
  sig_genes <- SKAT_results$SetID[which(SKAT_results$Q.value < 0.3)]
  
  ## subset sig genes from SSID to get table with gene & variants
  
  SSID_sig <- SSID[SSID$V1 %in% sig_genes,] 
  colnames(SSID_sig) <- c("gene", "variant")
  
  ## Calculate overall sample size
  N = length(phenotypes$phenotype)
  
  ## alpha = 0.05
  alpha = 0.05
  
  ## calculate power using two proportions (unequal n) test
  ##pwr.2p2n.test()
  
}

## Function to find matching values in a data frame and then return row and column names 
Which.names <- function(DF, value){
  ind <- which(DF==value, arr.ind=TRUE)
  row_name <- rownames(DF)[ind[,1]]
  col_name <- colnames(DF)[ind[,2]]
  row_N_cols <- data.frame(row_name, col_name)
  return(row_N_cols)
}

main()

