## Tiffany Timbers
## August 8, 2015
##
## Unblind measuremnet data from Volocity

main <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  data_file_input <- args[1]
  blind_key <- args[2]
  data_file_output <- args[3]
  
  stringsAsFactors = FALSE
  
  ## load libraries
  require(dplyr)
  require(stringr)
  
  ## read in data and key
  data <- read.table(data_file_input, header = FALSE, sep = ",", stringsAsFactors = FALSE)
  key <- read.table(blind_key, header = FALSE, sep = ",", stringsAsFactors = FALSE)
  
  ## make data header the header (Volocity csv's start with a blank line...)
  colnames(data) <- as.character(data[1,])
  data <- data[-1,]

  
  ## reduce data to necessary data and give meaningful column names
  measurement_ID <- data[1]
  key_ID <- data[2]
  length <- data[9]
  names(data)[1] <- "measurement_ID"
  names(data)[2] <- "key_ID"
  names(data)[9] <- "length"
  new_data <- data.frame(data$measurement_ID, data$key_ID, data$length)
  names(new_data)[1] <- "measurement_ID"
  names(new_data)[2] <- 'key_ID'
  names(new_data)[3] <- 'length'
  new_data$key_ID <- strtoi(new_data$key_ID)
  
  ## give column names to key
  colnames(key) <- c("key_ID", "filenames")
  
  ## make a new column called strain
  key$strain <- str_extract(key$filenames, "[A-Z]{2}[0-9]+")
  
  
  ## join new_data and key
  output_data <- left_join(new_data, key)

  write.table(output_data, data_file_output, row.names = FALSE, append = FALSE, quote = FALSE, sep = ',')
}

main()