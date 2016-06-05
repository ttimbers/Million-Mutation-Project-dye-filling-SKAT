# load libraries
library(stringr)
library(dplyr)
library(binom)
library(ggplot2)

# load data
data <- read.table("data/Crsipr_dye_filling.csv", sep = ',', header = TRUE, stringsAsFactors = FALSE)

# file prefix to save plot too
file_prefix1 <- "data/crispr_amphid"
file_prefix2 <- "data/crispr_phasmid"

# clean-up strain names
data$strain <- sub("3-12gr", "gk3639", data$strain)
data$strain <- sub("plasmid", "gk3637", data$strain)
data$strain <- sub("Cas9", "gk3674", data$strain)

# remove any spaces from "head" and "tail"
data$location <- sub("head ", "head", data$location)
data$location <- sub("tail ", "tail", data$location)

## Take data frame and aggregate over a list of the columns (e.g. strain, and time) and 
## returns a new data frame that can be used to plot probability with ggplot.

data_summary <- data %>% 
  group_by(location, strain) %>% 
  summarise(counts = sum(defect), N = length(defect), prob_dyf = sum(defect)/length(defect))

## Calculate the 95% binomial confidence intervals for the probabilities (Clopper-Pearson method)
conf_int <- binom.confint(data_summary$counts, data_summary$N, methods = "exact")

## Add these confidence intervals to the data frame
data_summary$conf_int_lower <- conf_int$lower
data_summary$conf_int_upper <- conf_int$upper

# make strain a factor
data_summary$strain <- as.factor(data_summary$strain)

# re-order factor levels
print(levels(data_summary$strain))
data_summary$strain = factor(data_summary$strain,levels(data_summary$strain)[c(5, 2, 3, 4, 1)])

## amphids
## make an object called my_plot which contains the plotting commands 
my_plot <- ggplot(data_summary[data_summary$location == "head",], aes(strain, prob_dyf)) + ## plot probability for each strain
  geom_bar(aes(group=strain), stat="identity") + ## make a line connecting all the points in the plot
  geom_errorbar(aes(ymin=conf_int_lower, ymax=conf_int_upper), ## add 95% confidence intervals
                width=.1) + ## make the confidence interval 0.1 width
  labs(x="Genotype", y="Proportion dye-filling \n defects") + ## label the x and y axes
  theme(plot.title = element_text(size = 16, vjust=2), ## Make the plot title larger and higher
        legend.title=element_blank(), ## remove the legend label
        legend.key=element_rect(fill='white'), ## remove the blocks around the legend items
        legend.text=element_text(size = 12), ## make the legend text font larger
        panel.background = element_rect(fill = 'grey96'), ## make the plot background grey
        axis.text.x=element_text(colour="black", size = 12, angle = 90), ## change the x-axis values   font to black and make larger
        axis.text.y=element_text(colour="black", size = 12), ## change the y-axis values font to black and make larger
        axis.title.x = element_text(size = 12, vjust = -0.2), ## change the x-axis label font to black, make larger, and move away from axis
        axis.title.y = element_text(size = 12, vjust = 1.3)) + ## change the y-axis label font to black, make larger, and move away from axis
  ylim(c(0,1)) ## Set the y-axis limits to a range from 0 to 1

## call the object to plot the figure
my_plot

ggsave(my_plot, filename = paste0(file_prefix1, '.pdf'), width = 3, height = 3.5, useDingbats=FALSE)

## phasmids
## make an object called my_plot which contains the plotting commands
my_plot <- ggplot(data_summary[data_summary$location == "tail",], aes(strain, prob_dyf)) + ## plot probability for each strain
  geom_bar(aes(group=strain), stat="identity") + ## make a line connecting all the points in the plot
  geom_errorbar(aes(ymin=conf_int_lower, ymax=conf_int_upper), ## add 95% confidence intervals
                width=.1) + ## make the confidence interval 0.1 width
  labs(x="Genotype", y="Proportion dye-filling \n defects") + ## label the x and y axes
  theme(plot.title = element_text(size = 16, vjust=2), ## Make the plot title larger and higher
        legend.title=element_blank(), ## remove the legend label
        legend.key=element_rect(fill='white'), ## remove the blocks around the legend items
        legend.text=element_text(size = 12), ## make the legend text font larger
        panel.background = element_rect(fill = 'grey96'), ## make the plot background grey
        axis.text.x=element_text(colour="black", size = 12, angle = 90), ## change the x-axis values   font to black and make larger
        axis.text.y=element_text(colour="black", size = 12), ## change the y-axis values font to black and make larger
        axis.title.x = element_text(size = 12, vjust = -0.2), ## change the x-axis label font to black, make larger, and move away from axis
        axis.title.y = element_text(size = 12, vjust = 1.3)) + ## change the y-axis label font to black, make larger, and move away from axis
  ylim(c(0,1)) ## Set the y-axis limits to a range from 0 to 1

## save the figure
ggsave(my_plot, filename = paste0(file_prefix2, '.pdf'), width = 3, height = 3.5, useDingbats=FALSE)

## Export stats results to a table
amphid_stats <- do_stats(data_summary, 'head')
write.table(amphid_stats, paste0(file_prefix1, '_stats.csv'), row.names = FALSE, append = FALSE, quote = FALSE, sep = ',')

phasmid_stats <- do_stats(data_summary, 'tail')
write.table(phasmid_stats, paste0(file_prefix2, '_stats.csv'), row.names = FALSE, append = FALSE, quote = FALSE, sep = ',')

# function to do stats
do_stats <- function(data_summary, column){
  ## make a dataframe containing the control strain phenotype count and N
  control.data  <- data_summary[data_summary$strain=="N2" & data_summary$location == column,]
  
  ## make wild-type/control contingency table
  defective <- control.data$counts
  wildtype <- control.data$N - control.data$counts
  control_data <- data.frame(cbind(defective, wildtype))
  
  ## remove control from phenotypes dataframe
  data_summary_muts <- data_summary[data_summary$strain != "N2" & data_summary$location == column,]
  
  ## Do a Fisher's to see if there is a difference between the groups
  ## make a dataframe to hold contingency table data for each strain
  defective <- data_summary_muts$counts
  wildtype <- data_summary_muts$N - data_summary_muts$counts
  strain <- as.character(data_summary_muts$strain_name)
  mutant_data <- data.frame(cbind(defective, wildtype))
  
  ## make a variable to hold p-values from the fisher's exact tests
  pvals <- rep(3,dim(data_summary_muts)[1])
  
  ## loop through and calculate p-value for each strain
  for (i in 1:dim(data_summary_muts)[1]) {
    row_1 <- control_data
    row_2 <- mutant_data[i,]
    contingency_table <- rbind(row_1, row_2) 
    temp <- fisher.test(contingency_table)
    index <- match(3,pvals)
    pvals[index] <- temp$p.value
  } 
  
  ## make a dataframe for p.values
  stat_results_table <- data.frame(data_summary_muts$strain, data_summary_muts$N, pvals)
  
  stat_results_table$test <- "fisher.test"
  stat_results_table$R_package <- "base R"
  
  ## add column of adjusted p-values
  stat_results_table$p_bonn_adjust <- p.adjust(stat_results_table$pvals, method = "bonferroni")
  
  return(stat_results_table)
}