## Tiffany Timbers
## 2015-07-28
## Analysis of Potential length defects

main <- function(){
  
  ## input_file <- "data/ADL_cilia_length_205-07-31.csv"
  
  ## Grab command line arguements
  args <- commandArgs(trailingOnly = TRUE)
  input_file <- args[1]
  what_was_measured <- args[2]
  figure_out <- args[3]
  stats_table_out <- args[4]
  
  ## load libraries
  require(ggplot2)
  require(stringr)
  require(plyr)
  require(binom)
  require(broom)

  stringsAsFactors = FALSE
  
  ## load data
  all_data <- read.csv(input_file, header=TRUE, sep=',', stringsAsFactors = FALSE)

  
  ## test for normality
  row1 <- tidy(shapiro.test(all_data$length[which(all_data$strain == "MX1924")]))
  row2 <- tidy(shapiro.test(all_data$length[which(all_data$strain == "MX2236")]))
  normality_stats <- rbind(row1, row2)
  stats_results <- normality_stats
  stats_results$test <- "shapiro.test"
  stats_results$Rpackage <- "base"
  stats_results$strain <- c("MX1924", "MX2236")
  
  
  ## test for difference using wilcox.test
  wilcox <- tidy(wilcox.test(all_data$length[which(all_data$strain == "MX1924")], all_data$length[which(all_data$strain == "MX2236")]))[1:2]
  wilcox$test <- "wilcox.test"
  wilcox$Rpackage <- "base"
  wilcox$strain <- c("MX1924 - MX2236")
  
  kruskal <- tidy(kruskal.test(all_data$length ~ as.factor(all_data$strain)))[1:2]
  kruskal$test <- "kruskal.test"
  kruskal$Rpackage <- "base"
  kruskal$strain <- c("MX1924 - MX2236")
  
  ## add to general stats data frame and save table
  stats_results <- rbind(stats_results, wilcox, kruskal)
  write.table(stats_results, stats_table_out, row.names = FALSE, append = FALSE, quote = FALSE)
  
  
  ## plot data
  g <- ggplot(all_data, aes(x = strain, y = length)) + ## plot lengths
    geom_boxplot(outlier.colour = NA) +
    theme(plot.title = element_text(size=20, face="bold", vjust=2), ## make the plot title larger and higher
          panel.background = element_rect(fill = "white"), ## make the plot background white
          axis.text.x=element_text(colour="black", size = 12), ## change the x-axis values font to black
          axis.text.y=element_text(colour="black", size = 12), ## change the y-axis values font to black and make larger
          axis.title.x = element_text(size = 12, vjust = -0.2), ## change the x-axis label font to black, make larger, and move away from axis
          axis.title.y = element_text(size = 12, vjust = 1.3)) +  ## change the y-axis label font to black, make larger, and move away from axis
    ##ggtitle("Violin Plot of Worm Area") +            ## set title
    labs(x="Genotype", y=paste(what_was_measured," length (um)", sep="")) +     ## label the x and y axes 
    geom_jitter(alpha = 0.7, position = position_jitter(width = 0.1), size = 2, colour="gray50")  ## overlay jitter plot
    
  #save plot
  ggsave(g, filename = figure_out, width = 4, height = 4)
  
  
}

main()