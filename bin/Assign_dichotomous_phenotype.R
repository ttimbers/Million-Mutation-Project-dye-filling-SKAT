## Do Fisher's exact test followed by a false-discovery rate (Benjamini-Hochberg 
## procedure) p-value adjustment to determine which strains are significantly different 
## from wild-type. 

main <- function(){
	args <- commandArgs(trailingOnly = TRUE)
	input_file_path <- args[1]
	output_file_path <- args[2]
	
	## read in data
	phenotype_counts <- read.csv(input_file_path, header=TRUE)
	
	## get dichotomous phenotypes for each strain from count data
	phenotypes.dichtomous <- get_dichotomous_phenotype(phenotype_counts[,1], phenotype_counts[,2], phenotype_counts[,3], "VC2010")
	
	## save this data to a file to be used for SKAT analysis
	write.table(phenotypes.dichtomous, output_file_path, sep="\t", row.names=FALSE, quote=FALSE, append=FALSE)
}

## Takes count data for a phenotype for a list of strains and uses a Fisher's exact test
## and and false-discovery rate (Benjamini-Hochberg procedure) p-value adjustment to 
## determine which strains are significantly different from wild-type. 
##
## Arguments:
## strain_name = a list of strain names
## defect_count = a list of the number animals with defects observed for each strain
## N = a list of thethe total number of animals assayed
## control_strain_name = a string that is the name of the strain to compare all others to
##
## Returns a dataframe of 3 columns, named "Strain" and  "Phenotype". Strain contains 
## strain names and "Phenotype" contains a 1 if strain's phenotype diverges significantly 
## from the control strain, and 0 if it does not. A 5% FDR (Benjamini-Hochberg procedure) 
## is used to adjust for multiple comparisons. 
##
get_dichotomous_phenotype  <- function(strain_name, defect_count, N, control_strain_name) {
  ## make a dataframe containing the strain_name, defect_count, and N
  counts <- data.frame(strain_name, defect_count, N)
  
  ## make a dataframe containing the control strain phenotype count and N
  control.data  <- counts[counts$strain_name==control_strain_name,]
  
  ## make wild-type/control contingency table
  defective <- control.data$defect_count
  wildtype <- control.data$N - control.data$defect_count
  control_data <- data.frame(cbind(defective, wildtype))

  rm(defective, wildtype)

  ## remove control from phenotypes dataframe
  counts <- subset(counts, strain_name != "VC2010")

  ## make a dataframe to hold contingency table data for each strain
  defective <- counts$defect_count
  wildtype <- counts$N - counts$defect_count
  strain <- as.character(counts$strain_name)
  mutant_data <- data.frame(cbind(defective, wildtype))

  ## make a variable to hold p-values from the fisher's exact tests
  pvals <- rep(3,dim(counts)[1])

  ## loop through and calculate p-value for each strain
  for (i in 1:dim(counts)[1]) {
    row_1 <- control_data
    row_2 <- mutant_data[i,]
    contingency_table <- rbind(row_1, row_2) 
    temp <- fisher.test(contingency_table)
    index <- match(3,pvals)
    pvals[index] <- temp$p.value
  } 

  ## make a dataframe for p.values
  fishers.results <- data.frame(counts$strain_name, pvals)

  ## FDR adjust p-values
  fdr.bh <- p.adjust(pvals, method = "BH", n = length(pvals))
  length(which(fdr.bh < 0.05))
  fdr.bh.results <- data.frame(strain, pvals, fdr.bh)

  ## classify which are fdr adjusted pvalues are < 0.05 as 1 (dyf) and those a > 0 as 0 (wt)
  fdr.bh.phenotypes <- fdr.bh.results
  fdr.bh.phenotypes$fdr.bh[fdr.bh.results$fdr.bh < 0.05] <- 1
  fdr.bh.phenotypes$fdr.bh[fdr.bh.results$fdr.bh > 0.05] <- 0

  fdr.bh.results$phenotype <- fdr.bh.phenotypes$fdr.bh

  colnames(fdr.bh.results)  <- c("strain", "pvals", "fdr.bh", "phenotype")
  return(fdr.bh.results)
}

main()
