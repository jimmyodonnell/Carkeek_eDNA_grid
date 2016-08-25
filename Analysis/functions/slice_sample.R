################################################################################
# A function to sample columns from a ( matrix | data frame | table )
# prop_prob : should the columns be sampled...
# ... randomly ('equal')
# ... proportional to their sum ('prop')
# ... inversely proportional to their sum ('invprop')
# frac : fraction of columns to sample
################################################################################

slice_sample <- function(dat, prop_prob = "equal", frac = 1) { # threshold
	
	if( class(dat) != "matrix" && class(dat) != "data.frame" ) {
		stop("'dat' must be a matrix or data frame")
	}
	
	if(prop_prob != "equal" && prop_prob != "prop" && prop_prob != "invprop" ) {
		stop("'prop_prob' must be 'equal', 'prop' or 'invprop'")		
	}
	
	if(frac > 1 | frac <= 0) {
		stop("'frac' must be > 0 and < = 1")
	}
	
	sums <- colSums(dat)
	
	if       (prop_prob == "equal") {
		probs <- NULL
	} else if(prop_prob == "prop") {
		probs <- sums/sum(sums)	
	} else if(prop_prob == "invprop" ) {
		probs <- (sum(sums)/sums)/sum(sum(sums)/sums)
	}

	samp_cols <- sample(
		x = 1:ncol(dat), 
		size = round(frac * ncol(dat)), 
		replace = FALSE, 
		prob = probs
		)
		
	return(dat[ , sort(samp_cols)])
	
}

# tests
# DAT <- matrix(data = rpois(100, lambda = 1), nrow = 10)
# slice_sample(dat = DAT, prop_prob = "equal", frac = 0.5)
