
no_turnover <- function(dataframe, grouping_vector)
{
	# this function calculates the maximum count
	# at which there is no longer turnover in presence/absence 
	# among subsamples (rows) within samples (grouping vector)
	
	# data should be a data frame before splitting
	the_data <- as.data.frame(dataframe)
	the_split <- split(the_data, f = grouping_vector)
	
	# for each subtable remove columns with sum = 0
	the_split <- lapply(the_split, strip_absent)
	
	# find the columns which are not > 0 in all rows
	temp <- lapply(the_split, not_in_all_rows)
	
	# calculate the maximum value
	return(sapply(temp, max))
	rm(temp)

}

# by_env_sample <- split(as.data.frame(otu_clean), metadata_clean[,colname_env_sample])
# by_env_sample <- lapply(by_env_sample, strip_absent)

# # find the columns which are not > 0 in all rows
# temp <- lapply(by_env_sample, not_in_all_rows)

# # calculate the maximum value
# sapply(temp, max)
# rm(temp)


# range(rowSums(otu_clean))
# range(rowSums(otu_scaled))