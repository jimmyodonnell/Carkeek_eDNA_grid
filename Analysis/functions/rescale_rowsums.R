rescale_rowsums <- function(mat)
{
	# calculate the minimum number of reads assigned to these OTUs in these samples
	minreads <- min(rowSums(mat))
	print(paste("minimum total counts per sample is", minreads))

	# scale the proportional abundance of the OTU in each sample to the minimum number of reads
	scaled <- prop(mat) * minreads

	# ignore counts OTUs found fewer than 0.5 times (anything greater would be rounded to 1)
	# otu_scaled[otu_scaled < 0.5] <- 0

	# round to whole numbers
	scaled <- round_right(scaled)
	dim(scaled)

	# exclude OTUs not found in these samples
	scaled <- strip_absent(scaled)
	
	# to re-order by the abundance in THESE samples (i.e. not those from samples from elsewhere)
	scaled <- scaled[,order(colSums(scaled), decreasing = TRUE)]

	print(paste("dimensions of rescaled matrix:", paste(dim(scaled), collapse = " ")))
	
	return(scaled)

}