library(qqman)

pvals <- read.table('data/Table_S4.csv', header = TRUE, sep = ' ')
qq(pvals$p.value)
pdf('data/qq_plot_amphid_log_plus_0_5.pdf')
qq(pvals$p.value, main = "Amphid dye-filling SKAT p-values \n using ln and a constant of 0.5")
dev.off()

pvals <- read.table('data/Table_S3.csv', header = TRUE, sep = ' ')
qq(pvals$p.value)
pdf('data/qq_plot_amphid_weights_log_plus_0_5.pdf')
qq(pvals$p.value, main = "Amphid dye-filling SKAT with weights p-values \n using ln and a constant of 0.5")
dev.off()

pvals <- read.table('/Users/tiffanytimbers/Documents/Post-Doc/Manuscripts/MMP_dyf_screen/PLoS_Genetics/tables/S4_Table.csv', header = TRUE, sep = ',')
qq(pvals$p.value)
pdf('data/qq_plot_amphid_logistic_regression.pdf')
qq(pvals$p.value, main = "Amphid dye-filling Logistic regression SKAT p-values")
dev.off()


pvals <- read.table('data/Table_S4.csv', header = TRUE, sep = ' ')
qq(pvals$p.value)
pdf('data/qq_plot_amphid_log10_plus_53-03.pdf')
qq(pvals$p.value, main = "Amphid dye-filling SKAT p-values \n using log10 and a constant of 5e-03")
dev.off()

pvals <- read.table('data/Table_S3.csv', header = TRUE, sep = ' ')
qq(pvals$p.value)
pdf('data/qq_plot_amphid_weights_log10_plus_53-03.pdf')
qq(pvals$p.value, main = "Amphid dye-filling SKAT with weights p-values \n using log10 and a constant of 5e-03")
dev.off()

pvals <- read.table('data/Table_S4.csv', header = TRUE, sep = ' ')
qq(pvals$p.value)
pdf('data/qq_plot_amphid_log10_plus_0_5.pdf')
qq(pvals$p.value, main = "Amphid dye-filling SKAT p-values \n using log10 and a constant of 0.5")
dev.off()

pvals <- read.table('data/Table_S3.csv', header = TRUE, sep = ' ')
qq(pvals$p.value)
pdf('data/qq_plot_amphid_weights_log10_plus_0_5.pdf')
qq(pvals$p.value, main = "Amphid dye-filling SKAT with weights p-values \n using log10 and a constant of 0.5")
dev.off()


# mystery analysis from May 15th
pvals <- read.table('/Users/tiffanytimbers/Downloads/Million-Mutation-Project-dye-filling-SKAT/data/Table_S4.csv', header = TRUE, sep = ' ')
qq(pvals$p.value)
pdf('data/qq_plot_amphid_log10_plus_0_5.pdf')
qq(pvals$p.value, main = "Amphid dye-filling SKAT p-values \n using log10 and a constant of 0.5")
dev.off()

pvals <- read.table('/Users/tiffanytimbers/Downloads/Million-Mutation-Project-dye-filling-SKAT/data/Table_S3.csv', header = TRUE, sep = ' ')
qq(pvals$p.value)
pdf('data/qq_plot_amphid_weights_log10_plus_0_5.pdf')
qq(pvals$p.value, main = "Amphid dye-filling SKAT with weights p-values \n using log10 and a constant of 0.5")
dev.off()
