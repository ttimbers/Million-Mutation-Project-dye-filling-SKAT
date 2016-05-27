# Script to plot power versus sample size from bootstrap power analysis

# Tiffany Timbers
# May 26, 2016

# Inputs:
#   1.  input_file - A tab separated file consisting of concatenated smaller files. It may include
#       repeated headers from the concatenation. It should have columns named SetID,	P.value,
#       N.Marker.All,	N.Marker.Test,	p_adjust,	N, ID
#
#   2.  output_file - name to save plot to
#
# Outputs:  A plot where the y-axis is power and the x-axis is sample size. Power is defined as the
#           proportion of bootstrap statistics obtained under the alternative hypothesis (a gene(s) is
#           signifcantly associated with the phenotype) that are more extreme than the critical values
#           (1 gene, 2 genes, 3 genes, etc.)
#

library(dplyr)
library(tidyr)
library(ggplot2)

input_file <- 'data/SKAT_power_summary_no_path.tsv'

# load the data
bootstrap_raw <- read.table(input_file, sep = '\t', header = TRUE)

# clean the data (may not always be necessary) to remove extra header rows
bootstrap_raw <- bootstrap_raw[bootstrap_raw$SetID != "SetID",]

# for each critical value, need to determine the power at each sample size
# map function to see if there is `critical_value` number of samples with p_adjust < 0.05
# if so, 1, if no, then 0

# nest/groupby/split by N
by_N <- bootstrap_raw %>% 
  group_by(N, ID) %>% 
  nest() %>% 
  mutate(one_gene = purrr::map(data, ~ significant(., 'p_adjust'))) %>% 
  mutate(two_genes = purrr::map(data, ~ significant(., 'p_adjust', 2))) %>% 
  mutate(three_genes = purrr::map(data, ~ significant(., 'p_adjust', 3))) %>% 
  mutate(four_genes = purrr::map(data, ~ significant(., 'p_adjust', 4))) %>% 
  mutate(five_genes = purrr::map(data, ~ significant(., 'p_adjust', 5))) %>% 
  unnest(one_gene) %>% 
  unnest(two_genes) %>% 
  unnest(three_genes) %>% 
  unnest(four_genes) %>% 
  unnest(five_genes)

# make long so it can be plotted
by_N_long <- gather(by_N, gene_n, sig, one_gene:five_genes)

# make gene_n a factor
by_N_long$gene_n <- as.factor(by_N_long$gene_n)

# get proportions to plot
power <- by_N_long %>% 
  group_by(N, gene_n) %>%  
  summarise(sum(sig)/length(sig))

colnames(power) <- c('N', 'gene_n', 'power')

# make N a number
power$N <- as.numeric(as.character(power$N))

# reorder levels

# make plot
power_plot <- ggplot(data = power, aes(x= N, y=power, group=gene_n, colour = gene_n)) +
  geom_line() +
  geom_point()


# save plot




# write a function that takes in a data frame and a column name, and returns a 1 if there are any 
# values less than alpha (default of alpha is 0.05) and a 0 if there are not
significant <- function(x, column, N = 1, alpha = 0.05) {
  # get values from rows where column is less than alpha
  sig_vals <- x[as.numeric(as.character(x[[column]])) < 0.05,]
  sig_vals <- sig_vals[[column]]
  if (length(sig_vals) >= N){
    return(1)
  } else {
    return(0)
  }
}
