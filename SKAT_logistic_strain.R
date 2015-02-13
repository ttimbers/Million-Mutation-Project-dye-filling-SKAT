##Load library for SKAT 
library(SKAT)

##Do SKAT on phasmid data (with and without weights)###############################################################
##
##Generate a SNP set data file (SSD) from binary plink formated data files using user specified SNP
##sets. 
Generate_SSD_SetID("phasmid/MMP_non-syn_coding.bed", "phasmid/MMP_non-syn_coding.bim", "phasmid/MMP_non-syn_coding.fam", "MMP_non-syn_coding_SSID.txt", "phasmid/MMP_phasmid.SSD", "phasmid/MMP_phasmid.info")

##read in fam file
fam_file <- read.table("phasmid/MMP_non-syn_coding.fam")

##put phenotypes into a vector
fam_phenotypes_vector <- fam_file$V6

##get mutation weights
SNPweights <- Read_SNP_WeightFile("MMP_SNP_WeightFile.txt")

##get SSD info from created file
SSD.info <- Open_SSD("phasmid/MMP_phasmid.SSD", "phasmid/MMP_phasmid.info")

##create null model based on phenotypes
Null_Model <- SKAT_Null_Model(fam_phenotypes_vector ~ 1, out_type="D")

##perform SKAT on all sets of variants (no weights)
All_SKAT_Data.no.weights  <- SKAT.SSD.All(SSD.info, Null_Model)
##sort All_SKAT_Data.no.weights by p-value
mydata.SKAT.no.weights <- All_SKAT_Data.no.weights$results
p.values.phasmid.no.weights <- mydata.SKAT.no.weights[order(mydata.SKAT.no.weights[,2]),]
write.table(p.values.phasmid.no.weights, "phasmid/SKAT_no_weights_phasmid_results", sep="\t", row.names=FALSE, quote=FALSE, append=FALSE)

##perform SKAT on all sets of variants with weights
All_SKAT_Data  <- SKAT.SSD.All(SSD.info, Null_Model, obj.SNPWeight=SNPweights)
##sort All_SKAT_Data by p-value
mydata.SKAT <- All_SKAT_Data$results
p.values.phasmid <- mydata.SKAT[order(mydata.SKAT[,2]),]
write.table(p.values.phasmid, "phasmid/SKAT_phasmid_weights_results", sep="\t", row.names=FALSE, quote=FALSE, append=FALSE)

##clean up environment
rm(list=ls())


##Do SKAT on amphid data (with and without weights)###############################################################
##
##Generate a SNP set data file (SSD) from binary plink formated data files using user specified SNP
##sets. 
Generate_SSD_SetID("amphid/MMP_non-syn_coding.bed", "amphid/MMP_non-syn_coding.bim", "amphid/MMP_non-syn_coding.fam", "MMP_non-syn_coding_SSID.txt", "amphid/MMP_amphid.SSD", "amphid/MMP_amphid.info")

##read in fam file
fam_file <- read.table("amphid/MMP_non-syn_coding.fam")

##put phenotypes into a vector
fam_phenotypes_vector <- fam_file$V6

##get mutation weights
SNPweights <- Read_SNP_WeightFile("MMP_SNP_WeightFile.txt")

##get SSD info from created file
SSD.info <- Open_SSD("amphid/MMP_amphid.SSD", "amphid/MMP_amphid.info")

##create null model based on phenotypes
Null_Model <- SKAT_Null_Model(fam_phenotypes_vector ~ 1, out_type="D")

##perform SKAT on all sets of variants (no weights)
All_SKAT_Data.no.weights  <- SKAT.SSD.All(SSD.info, Null_Model)
##sort All_SKAT_Data.no.weights by p-value
mydata.SKAT.no.weights <- All_SKAT_Data.no.weights$results
p.values.amphid.no.weights <- mydata.SKAT.no.weights[order(mydata.SKAT.no.weights[,2]),]
write.table(p.values.amphid.no.weights, "amphid/SKAT_no_weights_amphid_results", sep="\t", row.names=FALSE, quote=FALSE, append=FALSE)

##perform SKAT on all sets of variants with weights
All_SKAT_Data  <- SKAT.SSD.All(SSD.info, Null_Model, obj.SNPWeight=SNPweights)
##sort All_SKAT_Data by p-value
mydata.SKAT <- All_SKAT_Data$results
p.values.amphid <- mydata.SKAT[order(mydata.SKAT[,2]),]
write.table(p.values.amphid, "amphid/SKAT_amphid_weights_results", sep="\t", row.names=FALSE, quote=FALSE, append=FALSE)

##clean up environment
rm(list=ls())




##multiple testing corrections########################################################################################
##Adjust p-value for markers. Only considering markers with > 6 minor allele frequency (MAF)
library(fdrtool)

fdr_adjust_greater_6_markers <- function(phenotype.directory) {
  for (i in list.files(paste("./", phenotype.directory, sep=""),"*_results")) {
    df.p  <- read.table(paste("./", phenotype.directory, "/", i, sep=""), header=TRUE)
    df.p <- df.p[which(df.p$N.Marker.All > 6),]
    df.q <- fdrtool(df.p$P.value, statistic = "pvalue", cutoff.method="fndr")
    df.p$Q.value  <- df.q$qval
    write.table(df.p, paste("./", phenotype.directory, "/", i,"_w_qvalue", sep=""), sep="\t", row.names=FALSE, quote=FALSE, append=FALSE)
  }
}

fdr_adjust_greater_6_markers("amphid")
fdr_adjust_greater_6_markers("phasmid")