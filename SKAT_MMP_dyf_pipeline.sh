##Assign dichotomous phenotype to strains
Rscript Assign_dichotomous_phenotype.R

##make list_VCstrains_vcf.txt file for strains assayed
awk '{print $1}' dyf_phenotypes_dichotomous.txt | grep -h "^VC*" > vcf_files/list_VCstrains_vcf.txt

##Create one merged .vcf file for all strains that were assayed (listed in 
##a file called list_VCstrains_vcf.txt
##
##note: this .pl script also requires the use of a text file called 
##gknames.txt which contains all variant names of the MMP library, as 
##well as the VC*.combined.vcf files to combine into the merged file.
##These all should be in the same directory
cd vcf_files
./create_gwas_vcf_MMP.pl > ../MMP.vcf
cd ..

##create a merged .vcf file with only variants which affect coding regions
grep -v "GRANT=0" MMP.vcf > MMP_non-syn.vcf
grep -h "^#" MMP_non-syn.vcf > MMP_non-syn_coding.vcf
grep -h "coding" MMP_non-syn.vcf >> MMP_non-syn_coding_temp.vcf
grep -h -v "Noncoding" MMP_non-syn_coding_temp.vcf >> MMP_non-syn_coding.vcf
grep -h "GF=intron_splicing" MMP_non-syn.vcf >> MMP_non-syn_coding.vcf
rm MMP.vcf
rm MMP_non-syn.vcf

##make a version of the merged .vcf file without the header
Rscript Remove_vcf_header.R

##Make custom weights for each variant based off of the type of mutation
Rscript Make_custom_variant_weights.R

##Make snp set ID file (SSID) to run SKAT analysis
Rscript Make_SSID_file.R

##Make plink files from merged .vcf file to run SKAT analysis
Rscript Make_plink_files.R

##Do SKAT analysis
##Rscript SKAT_logistic_strain.R