################################################################################
# A function to remove from a ( matrix | data frame | table )
# columns which have a value in only one row.
################################################################################

rm_single_rows <- function(dat) {

	if(class(dat) != "matrix" && class(dat) != "data.frame") {
		stop("'dat' must be a matrix or data frame")
	}
	
	n_samples <- apply(X = dat, MARGIN = 2, FUN = function(x) sum(x > 0))
	
	return(dat[ , n_samples > 1])

}