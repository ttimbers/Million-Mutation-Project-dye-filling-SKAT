# Script to add a column of a repeated value to a dataframe
# Tiffany Timbers
# May 2016

# Inputs:
#   1. input file
#   2. value to add to column
#   3. column name
#   4. output file

main <- function(){
  
  # get command line args
  args <- commandArgs(trailingOnly = TRUE)
  input_file <- args[1]
  value <- args[2]
  column_name <- args[3]
  output_file <- args[4]
  
  # read in data
  data <- read.table(input_file, sep = '\t', header = TRUE)
  
  # add column of repeated value
  rep_value <- rep(x = value, dim(data)[1])
  data$new <- rep_value
  names(data)[names(data)=="new"] <- column_name
  
  # write file
  write.table(data, sep = '\t', col.names = TRUE, row.names = FALSE, append = FALSE, quote = FALSE)
}

main()