##Do Fisher's exact test to determine which strains are significantly different from wild-type. Wildtype has n=60, 0 amphid defects and 5 phasmid defects
wt.amphid <- data.frame(cbind(0,60))
colnames(wt.amphid) <- cbind("defective", "wt")
wt.phasmid <- data.frame(cbind(5,60))
colnames(wt.phasmid) <- cbind("defective", "wt")

##load phenotype data
phenotypes <- read.csv("/Users/michelleroux/Documents/Tiffany/SKAT/dyf_phenotpe_with_N.csv", header=TRUE)

##make a dataframe to hold contingency table data for each strain
amphid.defective <- round(phenotypes$Proportion.of.amphid.dye.fill.defects * phenotypes$N, digits=0)
amphid.wt <- phenotypes$N - amphid.defective
strain <- as.character(phenotypes$Strain)
phenotypes.amphid <- data.frame(cbind(strain, amphid.defective, amphid.wt))
colnames(phenotypes.amphid) <- cbind("strain", "defective", "wt")


phasmid.defective <- round(phenotypes$Proportion.of.phasmid.dye.fill.defects * phenotypes$N, digits=0)
phasmid.wt <- phenotypes$N - phasmid.defective
phenotypes.phasmid <- data.frame(cbind(strain, phasmid.defective, phasmid.wt))
colnames(phenotypes.phasmid) <- cbind("strain", "defective", "wt")

##make a variable to hold p-values from the fisher's exact tests
amphids.pvalues <- rep(3,480)
phasmids.pvalues <- rep(3,480)


##loop through and calculate p-value for each strain
for (i in 1:480) {
  amphid.row.1 <- cbind(wt.amphid[1,1], wt.amphid[1,2])
  amphid.row.2 <- cbind(as.numeric(levels(phenotypes.amphid[i,2]))[phenotypes.amphid[i,2]], as.numeric(levels(phenotypes.amphid[i,3]))[phenotypes.amphid[i,3]])
  amphid.contingency.table <- rbind(amphid.row.1, amphid.row.2) 
  amphid.x <- fisher.test(amphid.contingency.table)
  amphids.index <- match(3,amphids.pvalues)
  amphids.pvalues[amphids.index] <- amphid.x$p.value
  
  phasmid.row.1 <- cbind(wt.phasmid[1,1], wt.phasmid[1,2])
  phasmid.row.2 <- cbind(as.numeric(levels(phenotypes.phasmid[i,2]))[phenotypes.phasmid[i,2]], as.numeric(levels(phenotypes.phasmid[i,3]))[phenotypes.phasmid[i,3]])
  phasmid.contingency.table <- rbind(phasmid.row.1, phasmid.row.2) 
  phasmid.x <- fisher.test(phasmid.contingency.table)
  phasmids.index <- match(3,phasmids.pvalues)
  phasmids.pvalues[phasmids.index] <- phasmid.x$p.value
} 

##make a dataframe for p.values
fishers.results <- data.frame(strain,amphids.pvalues,phasmids.pvalues)

##how many strains will be significant if I use a conservative bonferroni correction:
length(which(fishers.results$amphids.pvalues < 0.0001))
length(which(fishers.results$phasmids.pvalues < 0.0001))

##FDR adjust p-values
fdr.bh.amphid <- p.adjust(amphids.pvalues, method = "BH", n = length(amphids.pvalues))
length(which(fdr.bh.amphid < 0.05))

fdr.bh.phasmid <- p.adjust(phasmids.pvalues, method = "BH", n = length(phasmids.pvalues))
length(which(fdr.bh.phasmid < 0.05))

fdr.bh.results <- data.frame(strain,fdr.bh.amphid,fdr.bh.phasmid)

##classify which are fdr adjusted pvalues are < 0.05 as 1 (dyf) and those a > 0 as 0 (wt)
fdr.bh.phenotypes <- fdr.bh.results
fdr.bh.phenotypes$fdr.bh.amphid[fdr.bh.results$fdr.bh.amphid < 0.05] <- 1
fdr.bh.phenotypes$fdr.bh.amphid[fdr.bh.results$fdr.bh.amphid > 0.05] <- 0
fdr.bh.phenotypes$fdr.bh.phasmid[fdr.bh.results$fdr.bh.phasmid < 0.05] <- 1
fdr.bh.phenotypes$fdr.bh.phasmid[fdr.bh.results$fdr.bh.phasmid > 0.05] <- 0

colnames(fdr.bh.phenotypes)  <- c("Strain", "dyf_amphid", "dyf_phasmid")
system("mkdir /Users/michelleroux/Documents/Tiffany/SKAT/strain_logistic")
write.table(fdr.bh.phenotypes, "/Users/michelleroux/Documents/Tiffany/SKAT/strain_logistic/dyf_phenotypes_dichotomous.txt", sep="\t", row.names=FALSE, quote=FALSE, append=FALSE)

rm(list=ls(all=TRUE))

##Make snp(variant) set ID (SSID) file for SKAT
##import vcf file without header
mydata.vcf <- read.table("/Users/michelleroux/Documents/Tiffany/SKAT/test_VCF_merge/MMPdyf_non-syn_coding_no_header.txt")

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
head(SSID)
tail(SSID)
SSID[10000:10010,]

##save SSID file
write.table(SSID, "/Users/michelleroux/Documents/Tiffany/SKAT/strain_logistic/MMPdyf_non-syn_coding_SSID.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, append=FALSE)

##remove everything but SSID
rm(variants, seq_tag, sequence_names, vars, seqs, mydata.vcf, SSID)

##create a folder for each phenotype
system("mkdir -p /Users/michelleroux/Documents/Tiffany/SKAT/strain_logistic/amphid")
system("mkdir -p /Users/michelleroux/Documents/Tiffany/SKAT/strain_logistic/phasmid")

##Run plink to make .ped and .fam files
system("/Users/michelleroux/Documents/plink_mac/plink --vcf /Users/michelleroux/Documents/Tiffany/SKAT/test_VCF_merge/MMPdyf_non-syn_coding.vcf --allow-extra-chr --out /Users/michelleroux/Documents/Tiffany/SKAT/strain_logistic/amphid/MMPdyf_non-syn_coding")
system("/Users/michelleroux/Documents/plink_mac/plink --vcf /Users/michelleroux/Documents/Tiffany/SKAT/test_VCF_merge/MMPdyf_non-syn_coding.vcf --allow-extra-chr --out /Users/michelleroux/Documents/Tiffany/SKAT/strain_logistic/phasmid/MMPdyf_non-syn_coding")

##add phenotypes to fam file
fam.file <- read.table("/Users/michelleroux/Documents/Tiffany/SKAT/strain_logistic/amphid/MMPdyf_non-syn_coding.fam")
fam.file <- fam.file[order(fam.file[,1]),]
phenotypes <- read.table("/Users/michelleroux/Documents/Tiffany/SKAT/strain_logistic/dyf_phenotypes_dichotomous.txt", header=TRUE)
phenotypes <- phenotypes[order(phenotypes[,1]),]
fam.file$V6  <- phenotypes$dyf_amphid
write.table(fam.file, "/Users/michelleroux/Documents/Tiffany/SKAT/strain_logistic/amphid/MMPdyf_non-syn_coding.fam", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, append=FALSE)
fam.file$V6  <- phenotypes$dyf_phasmid
write.table(fam.file, "/Users/michelleroux/Documents/Tiffany/SKAT/strain_logistic/phasmid/MMPdyf_non-syn_coding.fam", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, append=FALSE)
rm(fam.file, phenotypes, fam.ID)

rm(list=ls(all=TRUE))

##Make weights file for SKAT

##grep indels and snv's into two separate groups and create text files with this data
system("grep -h 'INDEL' /Users/michelleroux/Documents/Tiffany/SKAT/test_VCF_merge/MMPdyf_non-syn_coding_no_header.txt > /Users/michelleroux/Documents/Tiffany/SKAT/test_VCF_merge/MMPdyf_indels.txt")
system("grep -v 'INDEL' /Users/michelleroux/Documents/Tiffany/SKAT/test_VCF_merge/MMPdyf_non-syn_coding_no_header.txt > /Users/michelleroux/Documents/Tiffany/SKAT/test_VCF_merge/MMPdyf_snvs.txt")

##import text files created above
variants.indels <- read.table("/Users/michelleroux/Documents/Tiffany/SKAT/test_VCF_merge/MMPdyf_indels.txt")
variants.snvs <- read.table("/Users/michelleroux/Documents/Tiffany/SKAT/test_VCF_merge/MMPdyf_snvs.txt")


##Get variant names, amino acid change, indel label and reading frame change info for indels
##get indel variant names
variants.indels.names <- variants.indels$V3
vars.indel <- as.character(variants.indels.names)

##get indel amino acid change
library(stringr)
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

##Get variant names and reading frame change info for snvs
##get snv variant names
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
vars <- c(vars.indel, vars.snvs)
aacs  <- c(aacs.indel, aacs.snv)
indel.labels <- c(indel.label.indels, indel.label.snvs)
rfcs <- c(rfcs.indel, rfcs.snv)
coding.changes  <- data.frame(cbind(vars, aacs, indel.labels, rfcs))
head(coding.changes)
tail(coding.changes)
coding.changes[10000:10010,]

##save amino acid coding change file
write.table(coding.changes, "/Users/michelleroux/Documents/Tiffany/SKAT/strain_logistic/MMPdyf_coding_changes.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, append=FALSE)

##get inframe deletions (indel.lables == true and rfcs == no) and assign a weight of 0.5
indel.inframe <- subset(coding.changes, coding.changes$indel.labels == "yes" & coding.changes$rfcs == "no")
weight <- rep(x=0.5, dim(indel.inframe)[1])
indel.inframe <- data.frame(cbind(indel.inframe, weight))

##frameshift causing indels and assign a weight of 1
indel.frameshift <- subset(coding.changes, coding.changes$rfcs == "yes")
weight <- rep(x=1, dim(indel.frameshift)[1])
indel.frameshift <- data.frame(cbind(indel.frameshift, weight))

##nonsense mutations and assign a weight of 1
index.nonsense <- grep("[-][>][*]", coding.changes$aacs)
nonsense  <- coding.changes[index.nonsense,]
weight <- rep(x=1, dim(nonsense)[1])
nonsense <- data.frame(cbind(nonsense, weight))

##readthrough mutations and assign a weight of 1
index.readthrough <- grep("[*][-][>]", coding.changes$aacs)
readthrough  <- coding.changes[index.readthrough,]
weight <- rep(x=1, dim(readthrough)[1])
readthrough <- data.frame(cbind(readthrough, weight))

##indices of missense mutations and assign a weight of 0.25
index.muts <- grep("[-][>]", coding.changes$aacs)
muts <- coding.changes[index.muts,]
index.missense <- grep("[*]", muts$aacs, invert=TRUE)
missense  <- muts[index.missense,]
weight <- rep(x=0.25, dim(missense)[1])
missense <- data.frame(cbind(missense, weight))

##indices of splicing mutations and assign a weight of 1
splicing.defect <- subset(coding.changes, coding.changes$aacs == "NA")
weight <- rep(x=1, dim(splicing.defect)[1])
splicing.defect <- data.frame(cbind(splicing.defect, weight))

##combine all mutations
variant.weights <- rbind(splicing.defect, missense, readthrough, nonsense, indel.frameshift, indel.inframe)
vars.weights <- as.character(variant.weights$vars)
weight.weights <- as.character(variant.weights$weight)
variant.weights.file <- data.frame(cbind(vars.weights, weight.weights))

write.table(variant.weights.file, "/Users/michelleroux/Documents/Tiffany/SKAT/strain_logistic/MMP_SNP_WeightFile.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, append=FALSE)


##Load library for SKAT 
library(SKAT)

##Do SKAT on phasmid data with weights

##Generate a SNP set data file (SSD) from binary plink formated data files using user specified SNP
##sets. If you want to use plink formated data files, you must generate the SSD files first

##Make SSD for phasmid data
Generate_SSD_SetID("/Users/michelleroux/Documents/Tiffany/SKAT/strain_logistic/phasmid/MMPdyf_non-syn_coding.bed", "/Users/michelleroux/Documents/Tiffany/SKAT/strain_logistic/phasmid/MMPdyf_non-syn_coding.bim", "/Users/michelleroux/Documents/Tiffany/SKAT/strain_logistic/MMPdyf_non-syn_coding.fam", "/Users/michelleroux/Documents/Tiffany/SKAT/strain_logistic/MMPdyf_non-syn_coding_SSID.txt", "/Users/michelleroux/Documents/Tiffany/SKAT/strain_logistic/phasmid/MMPdyf_phasmid.SSD", "/Users/michelleroux/Documents/Tiffany/SKAT/strain_logistic/phasmid/MMPdyf_phasmid.info")

##read in fam file
##fam_file <- Read_Plink_FAM("/Users/michelleroux/Documents/Tiffany/SKAT/test_VCF_merge/phasmid/MMPdyf_phasmid.fam")
##above code gives me NA's in the phenotype column, and results in the ommission of 150 strains. Instead I will just use the read.table() function to open the fam_file
fam_file <- read.table("/Users/michelleroux/Documents/Tiffany/SKAT/strain_logistic/phasmid/MMPdyf_non-syn_coding.fam")

##put phenotypes into a vector
fam_phenotypes_vector <- fam_file$V6

##get mutation weights
SNPweights <- Read_SNP_WeightFile("/Users/michelleroux/Documents/Tiffany/SKAT/strain_logistic/MMP_SNP_WeightFile.txt")

##get SSD info from created file
SSD.info <- Open_SSD("/Users/michelleroux/Documents/Tiffany/SKAT/strain_logistic/phasmid/MMPdyf_phasmid.SSD", "/Users/michelleroux/Documents/Tiffany/SKAT/strain_logistic/phasmid/MMPdyf_phasmid.info")

##create null model based on phenotypes
Null_Model <- SKAT_Null_Model(fam_phenotypes_vector ~ 1, out_type="D")

##perform SKAT on all sets of variants (no weights)
All_SKAT_Data.no.weights  <- SKAT.SSD.All(SSD.info, Null_Model)
All_SKAT_Data.no.weights
str(All_SKAT_Data.no.weights)

##sort All_SKAT_Data by p-value
mydata.SKAT.no.weights <- All_SKAT_Data.no.weights$results
p.values.phasmid.no.weights <- mydata.SKAT.no.weights[order(mydata.SKAT.no.weights[,2]),]
write.table(p.values.phasmid.no.weights, "/Users/michelleroux/Documents/Tiffany/SKAT/strain_logistic/phasmid/SKAT_no_weights_phasmid_weights_results", sep="\t", row.names=FALSE, quote=FALSE, append=FALSE)

##perform SKAT on all sets of variants with weights
All_SKAT_Data  <- SKAT.SSD.All(SSD.info, Null_Model, obj.SNPWeight=SNPweights)
All_SKAT_Data
str(All_SKAT_Data)

##sort All_SKAT_Data by p-value
mydata.SKAT <- All_SKAT_Data$results
p.values.phasmid <- mydata.SKAT[order(mydata.SKAT[,2]),]
write.table(p.values.phasmid, "/Users/michelleroux/Documents/Tiffany/SKAT/strain_logistic/phasmid/SKAT_phasmid_weights_results", sep="\t", row.names=FALSE, quote=FALSE, append=FALSE)

##perform SKAT-O on all sets of variants
All_SKAT_Data.O  <- SKAT.SSD.All(SSD.info, Null_Model, obj.SNPWeight=SNPweights, method="optimal.adj")
All_SKAT_Data.O
str(All_SKAT_Data.O)

##sort All_SKAT_Data by p-value
mydata.SKAT.O <- All_SKAT_Data.O$results
p.values.phasmid.O <- mydata.SKAT.O[order(mydata.SKAT.O[,2]),]
write.table(p.values.phasmid.O, "/Users/michelleroux/Documents/Tiffany/SKAT/strain_logistic/phasmid/SKAT-O_phasmid_weights_results", sep="\t", row.names=FALSE, quote=FALSE, append=FALSE)

##do SKAT-O without weights
All_SKAT_Data.O.no.weights  <- SKAT.SSD.All(SSD.info, Null_Model, method="optimal.adj")
All_SKAT_Data.O.no.weights
str(All_SKAT_Data.O.no.weights)

##sort All_SKAT_Data.O.no.weights by p-value
mydata.SKAT.O.no.weights <- All_SKAT_Data.O.no.weights$results
p.values.phasmid.O.no.weights <- mydata.SKAT.O.no.weights[order(mydata.SKAT.O.no.weights[,2]),]
write.table(p.values.phasmid.O.no.weights, "/Users/michelleroux/Documents/Tiffany/SKAT/strain_logistic/phasmid/SKAT-O_no_weights_phasmid_weights_results", sep="\t", row.names=FALSE, quote=FALSE, append=FALSE)

##Do SKATBinary


##Do SKAT on amphid data with weights

##Generate a SNP set data file (SSD) from binary plink formated data files using user specified SNP
##sets. If you want to use plink formated data files, you must generate the SSD files first

##Make SSD for amphid data
Generate_SSD_SetID("/Users/michelleroux/Documents/Tiffany/SKAT/strain_logistic/amphid/MMPdyf_non-syn_coding.bed", "/Users/michelleroux/Documents/Tiffany/SKAT/strain_logistic/amphid/MMPdyf_non-syn_coding.bim", "/Users/michelleroux/Documents/Tiffany/SKAT/strain_logistic/amphid/MMPdyf_non-syn_coding.fam", "/Users/michelleroux/Documents/Tiffany/SKAT/strain_logistic/MMPdyf_non-syn_coding_SSID.txt", "/Users/michelleroux/Documents/Tiffany/SKAT/strain_logistic/amphid/MMPdyf_amphid.SSD", "/Users/michelleroux/Documents/Tiffany/SKAT/strain_logistic/amphid/MMPdyf_amphid.info")

##read in fam file
##fam_file <- Read_Plink_FAM("/Users/michelleroux/Documents/Tiffany/SKAT/test_VCF_merge/amphid/MMPdyf_amphid.fam")
##above code gives me NA's in the phenotype column, and results in the ommission of 150 strains. Instead I will just use the read.table() function to open the fam_file
fam_file.amphid <- read.table("/Users/michelleroux/Documents/Tiffany/SKAT/strain_logistic/amphid/MMPdyf_non-syn_coding.fam")

##put phenotypes into a vector
fam_phenotypes_vector.amphid <- fam_file.amphid$V6

##get mutation weights
SNPweights <- Read_SNP_WeightFile("/Users/michelleroux/Documents/Tiffany/SKAT/strain_logistic/MMP_SNP_WeightFile.txt")

##get SSD info from created file
SSD.info.amphid <- Open_SSD("/Users/michelleroux/Documents/Tiffany/SKAT/strain_logistic/amphid/MMPdyf_amphid.SSD", "/Users/michelleroux/Documents/Tiffany/SKAT/strain_logistic/amphid/MMPdyf_amphid.info")

##create null model based on phenotypes
Null_Model.amphid <- SKAT_Null_Model(fam_phenotypes_vector.amphid ~ 1, out_type="D")

##perform SKAT on all sets of variants
All_SKAT_Data.amphid  <- SKAT.SSD.All(SSD.info.amphid, Null_Model.amphid, obj.SNPWeight=SNPweights)
All_SKAT_Data.amphid
str(All_SKAT_Data.amphid)

##sort All_SKAT_Data by p-value
mydata.SKAT.amphid <- All_SKAT_Data.amphid$results
p.values.amphid <- mydata.SKAT.amphid[order(mydata.SKAT.amphid[,2]),]
write.table(p.values.amphid, "/Users/michelleroux/Documents/Tiffany/SKAT/strain_logistic/amphid/SKAT_amphid_weights_results", sep="\t", row.names=FALSE, quote=FALSE, append=FALSE)

##perform SKAT-O on all sets of variants
All_SKAT_Data.amphid.O  <- SKAT.SSD.All(SSD.info.amphid, Null_Model.amphid, obj.SNPWeight=SNPweights, method="optimal.adj")
All_SKAT_Data.amphid.O
str(All_SKAT_Data.amphid.O)

##sort All_SKAT_Data.O by p-value
mydata.SKAT.amphid.O <- All_SKAT_Data.amphid.O$results
p.values.amphid.O <- mydata.SKAT.amphid.O[order(mydata.SKAT.amphid.O[,2]),]
write.table(p.values.amphid.O, "/Users/michelleroux/Documents/Tiffany/SKAT/strain_logistic/amphid/SKAT-O_amphid_weights_results", sep="\t", row.names=FALSE, quote=FALSE, append=FALSE)





##work in progress... multiple testing corrections, display of p-values, calculating lambda (genomic inflation factor)

df.phasmid <- read.table("/Users/michelleroux/Documents/Tiffany/SKAT/strain_logistic/phasmid/SKAT_phasmid_weights_results", header=TRUE)

##FDR adjust p-values
p <- df.phasmid$P.value[which(df.phasmid$N.Marker.All > 3)]
fdr.SKAT.bh.phasmid <- p.adjust(p, method = "BH", n = length(p))
length(which(fdr.SKAT.bh.phasmid < 0.05))
min(fdr.SKAT.bh.phasmid)

library(qqman)
qq(df$P.value[which(df$N.Marker.All > 0)], main = "Q-Q plot of GWAS p-values")
qq(p.values.amphid$P.value[which(p.values.amphid$N.Marker.All > 1)], main = "Q-Q plot of GWAS p-values")

qqplot(ppoints(df$P.value[which(df$N.Marker.All > 2)]), df$P.value[which(df$N.Marker.All > 2)], main="Q-Q Plot of SKAT p-values", ylab="Sample Quantiles")

qqnorm(df$P.value[which(df$N.Marker.All > 2)])

hist(df$P.value[which(df$N.Marker.All > 1)], breaks = 10)


##calculate the genomic inflation factor, also known as lambda(λ)
lambdas <- rep(3, 30)
for (i in 1:30){
  set.seed(1121)
  pvalue1 <- df$P.value[which(df$N.Marker.All > i-1)]
  chisq1 <- qchisq(1-pvalue1,1)
  lambdas[i] = median(chisq1)/qchisq(0.5,1)
}

plot(lambdas)

##calculate the genomic inflation factor, also known as lambda(λ)
library(GenABEL)
lambdas.GenABEL <- rep(3, 30)
lambdas.GenABEL.CI  <- rep(3, 30)
for (i in 1:30){
    x <- estlambda(df$P.value[which(df$N.Marker.All > i-1)])
    lambdas.GenABEL[i]  <-  x$estimate
    lambdas.GenABEL.CI[i] <- x$se
}

df.lambda <- data.frame(lambdas.GenABEL, lambdas.GenABEL.CI)

library(ggplot2)
gp <- ggplot(df.lambda, aes(x=1:30, y=lambdas.GenABEL))
gp +  geom_point(size=2) + geom_errorbar(aes(ymax=lambdas.GenABEL+lambdas.GenABEL.CI, ymin=lambdas.GenABEL-lambdas.GenABEL.CI), width=.1) + xlab("Minimum minor allele count included") + ylab("Genomic inflation factor (lambda)") + ggtitle("Phasmid dye-filling SKAT with variant weights")

##calculate the genomic inflation factor, also known as lambda(λ)
set.seed(1234)
pvalue <- runif(1000, min=0, max=1)
data <- qchisq(pvalue, 1, lower.tail = FALSE)
data <- sort(df$P.value[which(df$N.Marker.All > 0)])
ppoi <- ppoints(data) #Generates the sequence of probability points
ppoi <- sort(qchisq(ppoi, df = 1, lower.tail = FALSE))
out <- list()
s <- summary(lm(data ~ 0 + ppoi))$coeff
out$estimate <- s[1, 1] # lambda 
out$se <- s[1, 2]
# median method
out$estimate <- median(data, na.rm = TRUE)/qchisq(0.5, 1)
