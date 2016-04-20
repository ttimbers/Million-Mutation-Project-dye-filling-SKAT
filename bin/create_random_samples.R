# Tiffany Timbers
# April 20, 2016
#
# Rscript that randomly samples, without replacement, samples (rows) from a csv file and 
# then saves the output of that file
#
# Arguments:
# 	1. input file
#   2. delimiter, can be either a comma (,) or a tab (/t)
#   3. TRUE/FALSE for whether or not there is a header
# 	4. column to randomly sample (number or name, and if name, header is assumed!)
#   5. number of samples to return
#	  6. output file
#
# Output:
#   1. a file containing the same columns as the input file, but with only N random samples

main <- function() {

	# get commmand line args
	args <- commandArgs(trailingOnly = TRUE)
	input_file <- args[1]
	delim <- args[2]
	print(delim)
	is_header <- as.logical(args[3])
	sample_col <- args[4]
	sample_N <- as.numeric(args[5])
	output_file	<- args[6]
	
	# load libraries
	library(stringr)
	
	# load input file
	if (str_detect(delim, 't')){
	  df <- read.table(file = input_file, header = is_header)
	} else if (str_detect(delim, ',')) {
	  df <- read.table(file = input_file, sep = ',', header = is_header)
	}
	print(head(df))
	
	# get random samples from the data frame
	sampled_df <- get_random_sample(x = df, col = sample_col, N = sample_N)
	print(head(sampled_df))
	
	# write sampled data frame to file
	if (str_detect(delim, 't')){
	  write.table(sampled_df, file = output_file, sep = '\t', col.names = is_header, row.names = FALSE, quote = FALSE, append = FALSE)
	} else if (str_detect(delim, ',')) {
	  write.table(sampled_df, file = output_file, sep = ',', col.names = is_header, row.names = FALSE, quote = FALSE, append = FALSE)
	}

}

# Randomly samples (without replacement) a data frame, and returns a smaller, randomly 
# sampled data frame. Assumes that you are indexing by column number if col begins with 
# the numbers between 1-9, otherwise, it assumes you are giving a column name.
# 
# Arguments:
#   1. x = a dataframe
#   2. col = column containing sample IDs
#   3. N = number of samples to extract
#
# Output
# a data frame with N random samples of x, from the column specified by col

get_random_sample <- function(x, col, N) {
  # sample data frame if number is provided for column
  if (str_detect(col, '[1-9]{1}')) {
    sampled_x <- x[sample(x[,as.numeric(col)], N, replace = FALSE, prob = NULL),]
    # sample data frame if name is provided for column    
  } else {
    sampled_x <- x[sample(x[,col], N, replace = FALSE, prob = NULL),]
  }

  return(sampled_x)  
}

main()