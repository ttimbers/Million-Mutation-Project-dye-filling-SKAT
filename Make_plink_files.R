##make plink files and add dichotomous phenotype data to plink fam file

make_plink <- function(vcf.filename, phenotypes_file, phenotype_column, where_to_save) {
  ##Takes a merged .vcf file and makes plink files from that .vcf file. It then takes a given phenotypes file and 
  ##alters the fam file to include the phenotypes
  ##
  ##Arguements:
  ##  vcf.filename: a .vcf file of all the strains you want to do SKAT on
  ##  phenotypes_file: a .txt file with a column listing strains, and subsequent columns listing phenotypes
  ##  phenotype_column: which column from phenotypes_file you want to grab the phenotypes from
  ##  where_to_save: a path to a directory to where you would like to output the plink files (including the altered fam file)
  ##
  ##Run plink to make .ped and .fam files
  merged_vcf_file_no_extension  <- substr(vcf.filename, 1, nchar(vcf.filename)-4)
  system(paste("plink --vcf ", vcf.filename," --allow-extra-chr --out ", where_to_save, "/", merged_vcf_file_no_extension, sep=""))
  
  ##add phenotypes to fam file
  fam.file <- read.table(paste(where_to_save, "/", merged_vcf_file_no_extension,".fam", sep=""))
  fam.file <- fam.file[order(fam.file[,1]),]
  phenotypes <- read.table(phenotypes_file, header=TRUE)
  phenotypes <- phenotypes[order(phenotypes[,1]),]
  fam.file$V6  <- phenotypes[,phenotype_column]
  write.table(fam.file, paste(where_to_save, "/", merged_vcf_file_no_extension,".fam", sep=""), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, append=FALSE)
}

##make directory to house amphid data and one to house phasmid data
system("mkdir amphid")
system("mkdir phasmid")

##make plink files and write amphid phenotype to fam file
make_plink("MMP_non-syn_coding.vcf", "dyf_phenotypes_dichotomous.txt", 2, "amphid")

vcf.filename  <- "MMP_non-syn_coding.vcf"
phenotypes_file <- "dyf_phenotypes_dichotomous.txt"
phenotype_column <- 2
where_to_save <- "amphid"

##make plink files and write amphid phenotype to fam file
make_plink("MMP_non-syn_coding.vcf", "dyf_phenotypes_dichotomous.txt", 3, "phasmid")