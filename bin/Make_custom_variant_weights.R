## Make weights file for SKAT
## Tiffany Timbers
## August 3, 2015
##
## Inputs:
##	1. .vcf file
## 	2. Weight to assign knockout mutations
## 	3. Weight to assign inframe deletion mutations
## 	4. Weight to assign missense mutations 
##	5. Path and filename of output file
##
## Outputs: a tab delimited text file of two columns. The first column is the variants
## and the second column is the weight assigned to the variant. 

main <- function() {
	
	## Start the clock!
  	ptm <- proc.time() 
	
	## load libraries
	require(stringr)

	## assign command line arguements
	args <- commandArgs(trailingOnly = TRUE)
	vcf.file <- args[1]
	weight.KO <- args[2]
	weight.inframe <- args[3]
  	weight.missense <- args[4]
	weights.file <- args[5]
	
	
	variant.weights <- custom_weights(vcf.file, weight.KO, weight.inframe, weight.missense)
	write.table(variant.weights, weights.file, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, append=FALSE)
	
	## output time to run script
  	the_time <- proc.time() - ptm # Stop the clock
  	print(paste("It took", the_time[3], "to run Make_custom_variant_weights.R")) 
}

custom_weights <- function(vcf.filename, weight.KO, weight.inframe.indel, weight.missense) {
  ## Arguments: 
  ## filename - name of a .vcf file with all of the strains you will be using for SKAT analysis
  ## weight.KO - custom weight to be assigned to knockout mutations (nonsense, splicing and frameshift causing indels)
  ## weight.inframe - custom weight to be assigned to indels that do not cause a frameshift
  ## weight.missense - custom weight to be assigned to missense mutations
  ##
  ## returns: a dataframe of two columns. The first column, vars.weight, is a list of variant names. The second
  ## column, weight.weights, is a list of corresponding weights assigned to the variants.

  #vcf_wout_header <- paste(substr(vcf.filename, 1, nchar(vcf.filename)-4), "_no_header.txt", sep="")
  vcf_beginning <- sub(".vcf", "", vcf.filename)
  
  ## grep indels and snv's into two separate groups and create text files with this data
  system(paste("grep -h 'INDEL' ", vcf.filename, " > ", vcf_beginning, ".indels", sep=""))
  system(paste("grep -h 'CODING=' ", vcf.filename, " > ", vcf_beginning, ".coding", sep=""))
  system(paste("grep -v 'INDEL' ", vcf.filename, " > ", vcf_beginning, ".snvs", sep="")) #GF=coding_exon

  ## import text files created above
  variants.indels <- read.table(paste(vcf_beginning, ".indels", sep=""))
  variants.coding <- read.table(paste(vcf_beginning, ".coding", sep=""))
  variants.snvs <- read.table(paste(vcf_beginning, ".snvs", sep=""))

  ## Get variant names, amino acid change, indel label and reading frame change info for indels
  ## get indel variant names
  variants.indels.names <- variants.indels$V3
  vars.indel <- as.character(variants.indels.names)

  ##get indel amino acid change
  aac_tag.indel  <- "AAC=[a-zA-Z0-9*>-_]{1,}"
  amino.acid.change.indel  <- str_extract(variants.indels$V8, aac_tag.indel)
  amino.acid.change.indel  <- sub("AAC=", "", amino.acid.change.indel)
  aacs.indel  <- as.character(amino.acid.change.indel)

  ##get indel label from indels
  indel.label.indels <- rep("yes", dim(variants.indels)[1])

  ##get reading frame change
  rfc_tag.indel  <- "RFC=[a-z]{1,}"
  reading.frame.change.indel  <- str_extract(variants.indels$V8, rfc_tag.indel)
  reading.frame.change.indel  <- sub("RFC=", "", reading.frame.change.indel)
  rfcs.indel <- as.character(reading.frame.change.indel)

  ## For .coding Get variant names, amino acid change, reading 
  ## frame change info
  variants.coding.names <- variants.coding$V3
  vars.coding <- as.character(variants.coding.names)
  
  ##get .coding amino acid change
  effect_tag.coding  <- "effect=[a-zA-Z_]{1,}"
  effect.coding  <- str_extract(variants.coding$V8, effect_tag.coding)
  effect.coding  <- sub("effect=", "", effect.coding)
  effect.coding  <- as.character(effect.coding)
  
  ##get indel label from coding
  indel.label.coding <- rep("yes", dim(variants.coding)[1])
  
  ##get reading frame change from coding
  rfcs.coding <- effect.coding
  
  ## Get variant names and reading frame change info for snvs
  ## get snv variant names
  variants.snvs.names <- variants.snvs$V3
  vars.snvs <- as.character(variants.snvs.names)

  ##get snv amino acid change
  aac_tag.snv  <- "AAC=[a-zA-Z0-9*>-]{1,}"
  amino.acid.change.snv  <- str_extract(variants.snvs$V8, aac_tag.snv)
  amino.acid.change.snv  <- sub("AAC=", "", amino.acid.change.snv)
  aacs.snv  <- as.character(amino.acid.change.snv)

  ##get indel label from snvs
  indel.label.snvs <- rep("no", dim(variants.snvs)[1])

  ##get reading frame change from snvs
  reading.frame.change.snv  <- rep("no", dim(variants.snvs)[1])
  rfcs.snv <- as.character(reading.frame.change.snv)

  ##concatenate indel variants with snv variants
  vars <- c(vars.indel, vars.coding, vars.snvs)
  aacs  <- c(aacs.indel, effect.coding, aacs.snv)
  indel.labels <- c(indel.label.indels, indel.label.coding, indel.label.snvs)
  rfcs <- c(rfcs.indel, rfcs.coding, rfcs.snv)
  coding.changes  <- data.frame(cbind(vars, aacs, indel.labels, rfcs))

  ##get inframe deletions (indel.lables == true and rfcs == no) and assign a specified weight
  indel.inframe <- subset(coding.changes, coding.changes$indel.labels == "yes" & coding.changes$rfcs == "no" | coding.changes$indel.labels == "yes" & coding.changes$rfcs == "noframeshift")
  weight <- rep(x=weight.inframe.indel, dim(indel.inframe)[1])
  indel.inframe <- data.frame(cbind(indel.inframe, weight))

  ##frameshift causing indels and assign a specified weight
  indel.frameshift <- subset(coding.changes, coding.changes$rfcs == "yes" | coding.changes$rfcs == "frameshift" | coding.changes$rfcs == "across_exon_boundary")
  weight <- rep(x=weight.KO, dim(indel.frameshift)[1])
  indel.frameshift <- data.frame(cbind(indel.frameshift, weight))

  ##nonsense mutations and assign a specified weight
  index.nonsense <- grep("[-][>][*]", coding.changes$aacs)
  nonsense  <- coding.changes[index.nonsense,]
  weight <- rep(x=weight.KO, dim(nonsense)[1])
  nonsense <- data.frame(cbind(nonsense, weight))

  ##readthrough mutations and assign a specified weight
  index.readthrough <- grep("[*][-][>]", coding.changes$aacs)
  readthrough  <- coding.changes[index.readthrough,]
  weight <- rep(x=weight.KO, dim(readthrough)[1])
  readthrough <- data.frame(cbind(readthrough, weight))

  ##indices of missense mutations and assign a specified weight
  index.muts <- grep("[-][>]", coding.changes$aacs)
  muts <- coding.changes[index.muts,]
  index.missense <- grep("[*]", muts$aacs, invert=TRUE)
  missense  <- muts[index.missense,]
  weight <- rep(x=weight.missense, dim(missense)[1])
  missense <- data.frame(cbind(missense, weight))

  ##indices of splicing mutations and assign a specified weight
  splicing.defect <- subset(coding.changes, coding.changes$aacs == "NA")
  weight <- rep(x=weight.KO, dim(splicing.defect)[1])
  splicing.defect <- data.frame(cbind(splicing.defect, weight))

  ##combine all mutations
  variant.weights <- rbind(splicing.defect, missense, readthrough, nonsense, indel.frameshift, indel.inframe)
  vars.weights <- as.character(variant.weights$vars)
  weight.weights <- as.character(variant.weights$weight)
  variant.weights.file <- data.frame(cbind(vars.weights, weight.weights))
  return(variant.weights.file)
}

main()