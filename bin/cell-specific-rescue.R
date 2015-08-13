## Tiffany Timbers
## 2015-08-12
## Analysis of cell-specific rescue of bgnt-1 cDNA in bgnt-1 mutants

main <- function(){
  
  ## input_file <- "data/cell-specific-rescue-amphid.csv"
  ## input_file <- "data/cell-specific-rescue-phasmid.csv"
  
  ## Grab command line arguements
  args <- commandArgs(trailingOnly = TRUE)
  input_file <- args[1]
  figure_out <- args[2]
  stats_table_out <- args[3]
  
  ## load libraries
  require(ggplot2)
  require(stringr)
  require(plyr)
  require(binom)
  require(broom)
  library(multcomp)
  
  ## load data
  all_data <- read.csv(input_file, header=TRUE)

  
  ## make data frame and column naming scheme more general
  col_to_agg_over <- c("strain")
  names(all_data)[names(all_data) == 'N_defects'] <- 'col_to_count'
  x <- all_data
  
  
  
  ## Take data frame and aggregate over a list of the columns (e.g. strain, and time) and 
  ## returns a new data frame that can be used to plot probability with ggplot.
  ##  
  ## New data frame contains columns: counts, N, prob, conf_int_lower and conf_int_upper. Confidnce 
  ## intervals are exact binomial confidence intervals (Pearson-Clopper method)
  ## For each strain, sum the col_to_count(s) and the sample size of that/those column(s)
  x_aggregate <- ddply(x, col_to_agg_over, summarise, counts = sum(col_to_count), N = sum(N))
  
  ## Calculate the probability for the col_to_count(s)
  x_aggregate <- ddply(x_aggregate, as.character(col_to_agg_over), transform, prob= (counts / N))
  
  ## Calculate the 95% binomial confidence intervals for the probabilities (Clopper-Pearson method)
  conf_int <- binom.confint(x_aggregate$counts, x_aggregate$N, methods = "exact")
  
  ## Add these confidence intervals to the data frame
  x_aggregate$conf_int_lower <- conf_int$lower
  x_aggregate$conf_int_upper <- conf_int$upper
  
  ## reorder levels of strain factor
  print(levels(x_aggregate$strain))
  x_aggregate$strain = factor(x_aggregate$strain,levels(x_aggregate$strain)[c(7, 6, 3, 4, 2, 1, 5)])
  
  
  
  ## make an object called my_plot which contains the plotting commands
  my_plot <- ggplot(x_aggregate, aes(strain, prob)) + ## plot probability for each strain
    geom_bar(aes(group=strain)) + ## make a line connecting all the points in the plot
    geom_point(size = 3) + ## make points larger
    geom_errorbar(aes(ymin=conf_int_lower, ymax=conf_int_upper), ## add 95% confidence intervals
                  width=.1) + ## make the confidence interval 0.1 width
    ##ggtitle('Habituation') + ## add a title to the plot
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
  
  ggsave(my_plot, filename = figure_out, width = 4, height = 4)
  
  ## make a dataframe containing the control strain phenotype count and N
  control.data  <- x_aggregate[x_aggregate$strain=="N2",]
  
  ## make wild-type/control contingency table
  defective <- control.data$counts
  wildtype <- control.data$N - control.data$counts
  control_data <- data.frame(cbind(defective, wildtype))
  
  ## remove control from phenotypes dataframe
  x_aggregate <- subset(x_aggregate, strain != "N2")
  
  ## Do a Fisher's to see if there is a difference between the groups
  ## make a dataframe to hold contingency table data for each strain
  defective <- x_aggregate$counts
  wildtype <- x_aggregate$N - x_aggregate$counts
  strain <- as.character(x_aggregate$strain_name)
  mutant_data <- data.frame(cbind(defective, wildtype))
  
  ## make a variable to hold p-values from the fisher's exact tests
  pvals <- rep(3,dim(x_aggregate)[1])
  
  ## loop through and calculate p-value for each strain
  for (i in 1:dim(x_aggregate)[1]) {
    row_1 <- control_data
    row_2 <- mutant_data[i,]
    contingency_table <- rbind(row_1, row_2) 
    temp <- fisher.test(contingency_table)
    index <- match(3,pvals)
    pvals[index] <- temp$p.value
  } 
  
  ## make a dataframe for p.values
  stat_results_table <- data.frame(x_aggregate$strain, pvals)
  stat_results_table$test <- "fisher.test"
  stat_results_table$R_package <- "base R"
  
  ## Export stats results to a table 
  write.table(stat_results_table, stats_table_out, row.names = FALSE, append = FALSE, quote = FALSE)
}

main()
