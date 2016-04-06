# dissimilarity by distance for reduced data set (mean abundance within sample for each OTU)

library(geosphere) # distm()

my_otu <- otu_mean
rownames(otu_mean) # PCT-C-0000 etc, aka "env_sample_name"
my_metadata <- metadata_exp[!duplicated(metadata_exp[,"env_sample_name"]),]

# calculate pairwise great circle distance between sampling locations using Haversine method
geo_dist <- distm(x = my_metadata[,c("lon", "lat")])
dimnames(geo_dist) <- list(my_metadata$env_sample_name, my_metadata$env_sample_name)

# test:
for(i in 1:length(my_metadata$env_sample_name)){
	if(geo_dist[my_metadata$env_sample_name[i], my_metadata$env_sample_name[i]] != 0){
		warning("distance between a point and itself is not zero, this indicates a problem")
	}
}

comm_dist <- as.matrix(vegdist(my_otu, method = "bray", diag = TRUE, upper = TRUE))
if(!(identical(dimnames(comm_dist), dimnames(geo_dist)))){
	warning("Whoa there! the row/column names of the two distance matrices do not seem to add up. This is bad.")
}

pdf(file = file.path(fig_dir, "diss_by_dist.pdf"), width = 8, height = 3)
	par(mar = c(4,5,1,1))
	plot(
		x = log(geo_dist + 100), 
		y = comm_dist, 
		ylim = c(0,1), 
		xaxt = "n", 
		pch = 1, 
		las = 1, 
		cex = 1, 
		col = rgb(0,0,0), #,alpha = 0.1 
		xlab = "Distance between samples (meters)", 
		ylab = "Bray-Curtis dissimilarity"
	)
	axis(side = 1, at = unique(log(metadata$dist_from_shore + 100)), labels = unique(metadata$dist_from_shore))
dev.off()
# abline(v = unique(log(metadata$dist_from_shore + 100)))





transects <- c("PCT-S", "PCT-C", "PCT-N")
diss_by_transect <- list()
for(transect in 1:length(transects)){
	diss_by_transect[[transect]] <- comm_dist[grep(transects[transect], rownames(comm_dist)), grep(transects[transect], colnames(comm_dist))]
}

dist_by_transect <- list()
for(transect in 1:length(transects)){
	dist_by_transect[[transect]] <- geo_dist[grep(transects[transect], rownames(geo_dist)), grep(transects[transect], colnames(geo_dist))]
}

mycols <- c("red", "green", "blue")
# mycols <- c("rgb(1,0,0)", "rgb(0,1,0)", "rgb(0,0,1)")

pdf(file = file.path(fig_dir, "diss_by_dist_by_transect.pdf"), width = 8, height = 3)
par(mar = c(4,5,1,1))
	plot(
		x = log(geo_dist + 100), 
		y = comm_dist, 
		ylim = c(0,1), 
		xaxt = "n", 
		pch = "", 
		las = 1, 
		cex = 1, 
		col = rgb(0,0,0), #,alpha = 0.1 
		xlab = "Distance between samples (meters)", 
		ylab = "Bray-Curtis dissimilarity"
	)
	axis(side = 1, at = unique(log(metadata$dist_from_shore + 100)), labels = unique(metadata$dist_from_shore))

for( i in 1:length(transects)){
	points(
		x = log(dist_by_transect[[i]] + 100), 
		y = diss_by_transect[[i]], 
		ylim = c(0,1), 
		xaxt = "n", 
		las = 1, 
		pch = 19, 
		cex = 1, 
		col = mycols[i], #,alpha = 0.1 
		xlab = "Distance between samples (meters)", 
		ylab = "Bray-Curtis dissimilarity"
	)
}
dev.off()
