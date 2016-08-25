# some functions for manipulating OTU tables.

scale01 <- function (x) {
# this function rescales a numeric vector to 0 and 1
  (x-min(x))/(max(x)-min(x))
}

as.binary <- function (a_matrix) {
# this function rescales a matrix to contain only 0 and 1
	bin_mat <- a_matrix
	bin_mat[bin_mat > 0] <- 1
	return(bin_mat)
}

strip_absent <- function(x) {
# this function removes any columns for which the sum is <= 0
	return(x[,which(colSums(x) > 0)])
}

present_in_all_rows <- function(x) {
# this function removes columns which are not > 0 in all rows
	return(x[ , colSums(x > 0) >= nrow(x)])
}

not_in_all_rows <- function(x) {
# this function removes columns which are > 0 in all rows
	return(x[ , colSums(x > 0) < nrow(x)])
}

prop <- function(x) {
# convert a matrix of counts to proportions of the samples
	return(x/rowSums(x))
}
