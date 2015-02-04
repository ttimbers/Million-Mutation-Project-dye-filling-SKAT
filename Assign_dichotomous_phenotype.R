##Do Fisher's exact test followed by a false-discovery rate (Benjamini-Hochberg procedure) p-value adjustment to determine which strains 
##are significantly different from wild-type. 

get_dichot_pheno  <- function(filename, control_strain) {
  ##load phenotype data
  phenotypes <- read.csv(filename, header=TRUE)

  #set wild-type/control strain
  control.data  <- phenotypes[phenotypes$Strain==control_strain,]
  #make wild-type/control contingency table
  wt.amphid <- data.frame(cbind(round(control.data$Proportion.of.amphid.dye.fill.defects*control.data$N,digits=0),control.data$N))
  colnames(wt.amphid) <- cbind("defective", "wt")
  wt.phasmid <- data.frame(cbind(round(control.data$Proportion.of.phasmid.dye.fill.defects*control.data$N,digits=0),control.data$N))
  colnames(wt.phasmid) <- cbind("defective", "wt")

  ##remove control from phenotypes dataframe
  phenotypes <- subset(phenotypes, Strain != "VC2010")

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
  amphids.pvalues <- rep(3,dim(phenotypes)[1])
  phasmids.pvalues <- rep(3,dim(phenotypes)[1])

  ##loop through and calculate p-value for each strain
  for (i in 1:dim(phenotypes)[1]) {
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
  return(fdr.bh.phenotypes)
}

phenotypes.dichtomous <- get_dichot_pheno("dyf_phenotpe_with_N_and_wt.csv", "VC2010")
write.table(phenotypes.dichtomous, "dyf_phenotypes_dichotomous.txt", sep="\t", row.names=FALSE, quote=FALSE, append=FALSE)

