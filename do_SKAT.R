main <- function() {
	args <- commandArgs(trailingOnly = TRUE)
	vcf.file <- args[1]
	SSID.file.name <- args[2]
	phenotype.list <- args[3]
	phenotype.file <- args[4]
	weights.file <- args[5]
	
	## make SSID file for my set of MMP strains
	## note - takes some time to read in MMPdyf_non-syn_coding_no_header.txt (derived from merged .vcf file)
	my.SSID <- make_SSID(vcf.file)
	## save SSID file
	write.table(my.SSID, SSID.file.name, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, append=FALSE)
	
	## make directory to house results of analysis for each phenotype, make plink files and then
	## write phenotype to the fam file
	phenos <- read.table(phenotype.list, header=FALSE)
	count = 1
	for p in phenos {
		system(paste("mkdir ", p, sep=""))
		count = count+1
		
		## make plink files and write phenotype to fam file
		make_plink(vcf.file, phenotype.file, count, p)
		
	}
	
	## Load library for SKAT 
	library(SKAT)

	## Do SKAT on phenotypes data (with and without weights)##############################
	for p in phenos {	
		## Navigate to phenotype directory
	 	system(paste("cd ", p, sep=""))
	 
	 	## Generate a SNP set data file (SSD) from binary plink formated data files using 
		## user specified SNP sets.
	 	Generate_SSD_SetID("*.bed", "*.bim", "*.fam", "../MMP_non-syn_coding_SSID.txt", paste("MMP_", p, ".SSD", sep=""), paste("MMP_", p, ".info", sep=""))

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
		p.values.phasmid.no.weights <- mydata.SKAT.no.weights[order(mydata.SKAT.no.weights[,2]),]
		write.table(p.values.phasmid.no.weights, paste(p_, "SKAT_no_weights_results.txt", sep=""), sep="\t", row.names=FALSE, quote=FALSE, append=FALSE)

		## perform SKAT on all sets of variants with weights
		All_SKAT_Data  <- SKAT.SSD.All(SSD.info, Null_Model, obj.SNPWeight=SNPweights)
		
		## sort All_SKAT_Data by p-value
		mydata.SKAT <- All_SKAT_Data$results
		p.values.phasmid <- mydata.SKAT[order(mydata.SKAT[,2]),]
		write.table(p.values.phasmid, paste(p_, "SKAT_weights_results.txt", sep=""), sep="\t", row.names=FALSE, quote=FALSE, append=FALSE)

		##navigate to parent directory (i.e. out of phenotype directory
		system("cd ..")
		
		## do False Discovery Rate analysis
		fdr_adjust_greater_6_markers(p)
	}
}

## Make snp(variant) set ID (SSID) file for SKAT
make_SSID <- function(vcf.filename) {
  ## takes merged .vcf file and creates a dataframe from which you can save a snp(variant) set ID (SSID) 
  ## file for SKAT analysis
  ##
  vcf_sans_header <- paste(substr(vcf.filename, 1, nchar(vcf.filename)-4), "_no_header.txt", sep="")
  
  ## import vcf file without header
  mydata.vcf <- read.table(vcf_sans_header)

  ## get variant names
  variants <- mydata.vcf$V3

  ## get sequence names
  library(stringr)
  seq_tag  <- "SN=[a-zA-Z0-9.]{1,}"
  sequence_names  <- str_extract(mydata.vcf$V8, seq_tag)
  sequence_names  <- sub("SN=", "", sequence_names)

  vars <- as.character(variants)
  seqs  <- as.character(sequence_names)

  SSID  <- data.frame(cbind(seqs, vars))
  return(SSID)
}

## make plink files and add dichotomous phenotype data to plink fam file
make_plink <- function(vcf.filename, phenotypes_file, phenotype_column, where_to_save) {
  ## Takes a merged .vcf file and makes plink files from that .vcf file. It then takes a given phenotypes file and 
  ## alters the fam file to include the phenotypes
  ##
  ## Arguements:
  ##   vcf.filename: a .vcf file of all the strains you want to do SKAT on
  ##   phenotypes_file: a .txt file with a column listing strains, and subsequent columns listing phenotypes
  ##   phenotype_column: which column from phenotypes_file you want to grab the phenotypes from
  ##   where_to_save: a path to a directory to where you would like to output the plink files (including the altered fam file)
  ##
  ## Run plink to make .ped and .fam files
  merged_vcf_file_no_extension  <- substr(vcf.filename, 1, nchar(vcf.filename)-4)
  system(paste("plink --vcf ", vcf.filename," --allow-extra-chr --out ", where_to_save, "/", merged_vcf_file_no_extension, sep=""))
  
  ## add phenotypes to fam file
  fam.file <- read.table(paste(where_to_save, "/", merged_vcf_file_no_extension,".fam", sep=""))
  fam.file <- fam.file[order(fam.file[,1]),]
  phenotypes <- read.table(phenotypes_file, header=TRUE)
  phenotypes <- phenotypes[order(phenotypes[,1]),]
  fam.file$V6  <- phenotypes[,phenotype_column]
  write.table(fam.file, paste(where_to_save, "/", merged_vcf_file_no_extension,".fam", sep=""), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, append=FALSE)
}

## generate q-values using false discovery rate using fdrtool()
fdr_adjust_greater_6_markers <- function(phenotype.directory) {
  library(fdrtool)
  for (i in list.files(paste("./", phenotype.directory, sep=""),"*_results")) {
    df.p  <- read.table(paste("./", phenotype.directory, "/", i, sep=""), header=TRUE)
    df.p <- df.p[which(df.p$N.Marker.All > 6),]
    df.q <- fdrtool(df.p$P.value, statistic = "pvalue", cutoff.method="fndr")
    df.p$Q.value  <- df.q$qval
    write.table(df.p, paste("./", phenotype.directory, "/", i,"_w_qvalue", sep=""), sep="\t", row.names=FALSE, quote=FALSE, append=FALSE)
  }
}

main()