## Tiffany Timbers
## 2015-07-28
## Analysis of Potential guidance defects in bgnt-1 mutants

main <- function(){
  
  ## input_file <- "data/ADL_guidance.csv"
  
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

  ## load data
  all_data <- read.csv(input_file, header=TRUE)

  
  ## clean data
  ## remove spaces from strain names
  unique(all_data$strain)
  all_data$strain <- sub("MX1924 ", "MX1924", all_data$strain)
  
  ## make data frame just for data concerning > 1 double rod cilia
  mult_double_rod <- all_data[,1:6]
  mult_double_rod <- cbind(mult_double_rod, all_data[,8:11])

  ## remove NaNs
  mult_double_rod <- mult_double_rod[complete.cases(mult_double_rod), ]
  
  ## make data frame and column naming scheme more general
  col_to_agg_over <- c("strain")
  names(mult_double_rod)[names(mult_double_rod) == 'multiple_double_rod_cilia_on_one_side'] <- 'col_to_count'
  x <- mult_double_rod
  
  
  
  ## Take data frame and aggregate over a list of the columns (e.g. strain, and time) and 
  ## returns a new data frame that can be used to plot probability with ggplot.
  ##  
  ## New data frame contains columns: counts, N, prob, conf_int_lower and conf_int_upper. Confidnce 
  ## intervals are exact binomial confidence intervals (Pearson-Clopper method)
  ## For each strain, sum the col_to_count(s) and the sample size of that/those column(s)
  x_aggregate <- ddply(x, col_to_agg_over, summarise, counts = sum(col_to_count), N = length(col_to_count))
  
  ## Calculate the probability for the col_to_count(s)
  x_aggregate <- ddply(x_aggregate, as.character(col_to_agg_over), transform, prob= (counts / N))
  
  ## Calculate the 95% binomial confidence intervals for the probabilities (Clopper-Pearson method)
  conf_int <- binom.confint(x_aggregate$counts, x_aggregate$N, methods = "exact")
  
  ## Add these confidence intervals to the data frame
  x_aggregate$conf_int_lower <- conf_int$lower
  x_aggregate$conf_int_upper <- conf_int$upper
  
  
  
  ## make an object called my_plot which contains the plotting commands
  my_plot <- ggplot(x_aggregate, aes(strain, prob)) + ## plot probability for each strain
  geom_bar(aes(group=strain)) + ## make a line connecting all the points in the plot
  geom_point(size = 3) + ## make points larger
  geom_errorbar(aes(ymin=conf_int_lower, ymax=conf_int_upper), ## add 95% confidence intervals
                width=.1) + ## make the confidence interval 0.1 width
  ##ggtitle('Habituation') + ## add a title to the plot
  labs(x="Genotype", y="Proportion with > 1 \n double rod cilia/amphid") + ## label the x and y axes
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
  
  
  
  ## Do a fischer's exact test to see if there is a difference between the groups
  col_1 <- x_aggregate[,2]
  col_2 <- x_aggregate[,3] - x_aggregate[,2]
  contingency_table <- cbind(col_1, col_2) 
  stat_results <- fisher.test(contingency_table)
  stat_results_table <- tidy(stat_results)
  stat_results_table$test <- "fisher.test"
  stat_results_table$R_package <- "base R"
  
  ## Export stats results to a table 
  write.table(stat_results_table, stats_table_out, row.names = FALSE, append = FALSE, quote = FALSE)
}

## Takes in a data frame, a list of the columns to aggregate over (e.g. strain, and time) and a list 
## of columns to count (and get probability for). This returns a new data frame that can be used  
## to plot probability with ggplot. 
## 
## New data frame contains columns: counts, N, prob, conf_int_lower and conf_int_upper. Confidnce 
## intervals are exact binomial confidence intervals (Pearson-Clopper method)
##
## Dependencies: plyr and binom
create_prob_table_to_plot <- function(x, col_to_agg_over, col_to_count){
  ## For each strain, sum the col_to_count(s) and the sample size of that/those column(s)
  x_aggregate <- ddply(x, as.character(col_to_agg_over), summarise, counts = sum(col_to_count), N = length(col_to_count))
  
  ## Calculate the probability for the col_to_count(s)
  x_aggregate <- ddply(x_aggregate, as.character(col_to_agg_over), transform, prob= (counts / N))
  
  ## Calculate the 95% binomial confidence intervals for the probabilities (Clopper-Pearson method)
  conf_int <- binom.confint(x_aggregate$counts, x_aggregate$N, methods = "exact")
  
  ## Add these confidence intervals to the data frame
  x_aggregate$conf_int_lower <- conf_int$lower
  x_aggregate$conf_int_upper <- conf_int$upper
  
  return(x_aggregate)
}

main()
  