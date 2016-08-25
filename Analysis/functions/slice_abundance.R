################################################################################
# A function to remove columns from a ( matrix | data frame | table )
# if the column sum is above/below a given threshold
################################################################################

slice_abundance <- function(dat, threshold, keep = "gt") { # threshold
	
	if(class(threshold) != "numeric") {
		stop("'threshold' must be numeric")
	}
	
	if(length(threshold) > 1) {
		warning("'threshold' has more than 1 value; only using the first")
		threshold <- threshold[1]
	}
	
	if(keep != "gt" && keep != "lt") {
		stop("'keep' must be either 'lt' or 'gt'")		
	}
	
	if(class(dat) != "matrix" && class(dat) != "data.frame") {
		stop("'dat' must be a matrix or data frame")
	}
	
	sums <- colSums(dat)
	
	if        (keep == "gt") {
		sliced <- dat[ , sums > threshold ]	
	} else if (keep == "lt") {
		sliced <- dat[ , sums < threshold ]	
	}
	
	return(sliced)
}

# tests
# slice_abundance(threshold = c(1,2))
# slice_abundance(threshold = "a")
# slice_abundance(threshold = 1, keep = "pt")

# DAT <- matrix(data = rpois(100, lambda = 1), nrow = 10)
# slice_abundance(dat = "cat", threshold = 2, keep = "gt")
# slice_abundance(dat = DAT, threshold = 2, keep = "gt")
# slice_abundance(dat = DAT, threshold = 8, keep = "gt")
# slice_abundance(dat = DAT, threshold = 10, keep = "lt")
