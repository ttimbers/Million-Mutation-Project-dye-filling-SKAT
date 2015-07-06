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

main <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  path_to_phenotypes <- args[1]
  path_to_SKAT_analysis <- args[2]
  path_to_SSID <- args[3]
  path_to_vcf_w_oneline_header <- args[4]
  
  require(pwr)
  require(dplyr)
  require(plyr)
  
  ## Read in phenotype file
  phenotypes <- read.csv(path_to_phenotypes, header = TRUE, sep = "\t")
  
  ## Read in SKAT analysis results file
  SKAT_results <- read.table(path_to_SKAT_analysis, header = TRUE, sep="\t")
  
  ## Read in SSID file
  SSID <- read.table(path_to_SSID, sep="\t")
  colnames(SSID) <- c("gene", "variant")
  
  ## Reduce SSID to only variants which were included in SKAT analysis
  genes_to_incl <- SKAT_results$SetID
  SSID <- SSID[SSID$gene %in% genes_to_incl,]
  
  ## Read in VCF file with oneline header
  VCF <- read.table(path_to_vcf_w_oneline_header, header=TRUE)
  
  ## Initialize variables to append to in loop
  strain <- c()
  add_on_after_strain <- c()
  add_on_after_variant <- c()
  add_on_after_gene <- c()
   
  ## Loop through variants and input strain name into SSID_sig
  for (i in 1:dim(SSID)[1]) {
  #for (i in 36) {
    var_row <- VCF[which(VCF$ID == SSID$variant[i]), ]
    var_strain <- names(which(apply(var_row, 2, function(x) any(grepl("1/1", x)))))
    #print(i)
    #print(var_strain)
    
    if (length(var_strain) == 2){
      add_on_after_strain <- c(add_on_after_strain, var_strain[2])
      add_on_after_variant <- c(add_on_after_variant, as.character(SSID$variant[i]))
      add_on_after_gene <- c(add_on_after_gene, as.character(SSID$gene[i]))} 
      strain <- c(strain, var_strain[1])
  }
  
  ## Add strain to SSID_sig
  SSID$strain <- strain
  
  ## Make a dataframe of strains that had duplicated variants
  add_on_after_all <- data.frame(add_on_after_strain, add_on_after_variant, add_on_after_gene) 
  colnames(add_on_after_all) <- c("gene", "variant", "strain")
  
  ## Concatenate these two data frames together
  SSID <- rbind(SSID, add_on_after_all)
  
  ## Add phenotypes to SSID
  phenotypes_reduced <- phenotypes[,1:2]
  df_for_pwr <- right_join(SSID, phenotypes_reduced, by = "strain")
  ##write.table(df_for_pwr, file=)
  
  ## Determine average effect size for significantly variants from SKAT analysis
  ##
  effect_size <- ddply(df_for_pwr, c("gene"), summarise,
                 N_variants  = length(variant),
                 N_cases = sum(phenotype))
  
  effect_size$N_controls <- effect_size$N_variants - effect_size$N_cases
  
  
 
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

main()

