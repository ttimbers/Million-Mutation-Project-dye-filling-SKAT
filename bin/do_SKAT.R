main <- function() {
	args <- commandArgs(trailingOnly = TRUE)
	path_to_fam_file <- args[1]
	path_to_phenotype.filename <- args[2]
	path_to_plink_files <- args[3]
  	SSID_file <- args[4]
	weights.file <- args[5]
	
	require(SKAT)
	require(fdrtool)
	
	## write phenotype to the fam file
	write_to_fam(path_to_fam_file, path_to_phenotype.filename)	

	## Do SKAT on phenotypes data (with and without weights)
	 
	## Generate a SNP set data file (SSD) from binary plink formated data files using 
	## user specified SNP sets.
  
	Generate_SSD_SetID(paste(path_to_plink_files, "/", list.files(path_to_plink_files, pattern='*.bed'), sep=""), paste(path_to_plink_files, "/", list.files(path_to_plink_files, pattern='*.bim'), sep=""), paste(path_to_plink_files, "/", list.files(path_to_plink_files, pattern='*.fam'), sep=""), SSID_file, paste(path_to_plink_files, "/file.SSD", sep=""), paste(path_to_plink_files, "/data/file.info", sep=""))
	
	## read in fam file
	fam_file <- read.table("*.fam")
		
	## put phenotypes into a vector
	fam_phenotypes_vector <- fam_file$V6
		
	## get mutation weights
	SNPweights <- Read_SNP_WeightFile(paste("../", weights.file, sep=""))
		
	## get SSD info from created file
	SSD.info <- Open_SSD("*.SSD", "*.info")
		
	## create null model based on phenotypes
	Null_Model <- SKAT_Null_Model(fam_phenotypes_vector ~ 1, out_type="D")

	## perform SKAT on all sets of variants (no weights)
	All_SKAT_Data.no.weights  <- SKAT.SSD.All(SSD.info, Null_Model)
		
	## sort All_SKAT_Data.no.weights by p-value
	mydata.SKAT.no.weights <- All_SKAT_Data.no.weights$results
	p.values.no.weights <- mydata.SKAT.no.weights[order(mydata.SKAT.no.weights[,2]),]
	write.table(p.values.no.weights, "SKAT_no_weights_results.txt", sep="\t", row.names=FALSE, quote=FALSE, append=FALSE)

	## do False Discovery Rate analysis for SKAT without weights
	pq_wo_weights <- fdr_adjust(p.values.no.weights, 6)
	write.table(pq_wo_weights, "SKAT_pANDq_no_weights_results.txt", sep="\t", row.names=FALSE, quote=FALSE, append=FALSE)

	## perform SKAT on all sets of variants with weights
	All_SKAT_Data  <- SKAT.SSD.All(SSD.info, Null_Model, obj.SNPWeight=SNPweights)
		
	## sort All_SKAT_Data by p-value
	mydata.SKAT.weights <- All_SKAT_Data$results
	p.values.weights <- mydata.SKAT.weights[order(mydata.SKAT.weights[,2]),]
	write.table(p.values.weights, "SKAT_weights_results.txt", sep="\t", row.names=FALSE, quote=FALSE, append=FALSE)

	## do False Discovery Rate analysis for SKAT with weights
	pq_weights <- fdr_adjust(p.values.weights, 6)
	write.table(pq_weights, "SKAT_pANDq_weights_results.txt", sep="\t", row.names=FALSE, quote=FALSE, append=FALSE)
	
	##navigate to parent directory (i.e. out of phenotype directory
	system("cd ..")
}

## Adds dichotomous phenotype data to plink fam file
write_to_fam <- function(path, phenotypes) {
  ##
  ## Arguements:
  ##   path: path to the plink fam file to be added to
  ##   phenotypes: a file where the first column lists strains and the second column lists 
  ##   phenotype
  ## 
  fam.file <- read.table(path)
  fam.file <- fam.file[order(fam.file[,1]),]
  pheno_vals <- read.table(phenotypes, header=TRUE)
  pheno_vals <- pheno_vals[order(pheno_vals[,1]),]
  fam.file$V6  <- pheno_vals[,2]
  write.table(fam.file, path, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, append=FALSE)
}

## generate q-values using false discovery rate using fdrtool()
fdr_adjust_greater_6_markers <- function(pvals, minNmarkers) {
    pvals <- pvals[which(pvals$N.Marker.All > minNmarkers),]
    qvals <- fdrtool(pvals$P.value, statistic = "pvalue", cutoff.method="fndr")
    pvals$Q.value  <- qvals$qval
    pAndqvals <- pvals
    return(pAndqvals) 
}

main()