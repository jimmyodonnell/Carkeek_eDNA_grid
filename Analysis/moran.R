# Calculate Moran's I

my_table    <- otu_mean[,1:50]
my_metadata <- metadata_mean

library(ape) # Moran.I

my_dist <-  distm(my_metadata[,c(colname_lon, colname_lat)])
my_dist_inv <- 1/my_dist
diag(my_dist_inv) <- 0

# Calculate Moran's I for the taxa:
moran_out <- t(
	apply(
		X = my_table, 
		MARGIN = 2, 
		FUN = function(x) {
			unlist(Moran.I(x = x, weight = my_dist_inv))
		}
		)
	)

plot(sort(moran_out[,"p.value"]))
abline(h = 0.05, lty = "dashed")

# is significance correlated with abundance?
plot(colSums(my_table), moran_out[,"p.value"])
abline(h = 0.05, lty = "dashed")

# is significance correlated with the number of sites at which a taxon occurs?
plot(colSums(my_table > 1), moran_out[,"p.value"])
abline(h = 0.05, lty = "dashed")

