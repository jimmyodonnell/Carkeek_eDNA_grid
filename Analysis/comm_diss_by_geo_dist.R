# FIRST RUN 'load_data.r'

library(geosphere) # distm()

# calculate pairwise great circle distance between sampling locations using Haversine method
geo_dist <- distm(x = metadata[,c("lon", "lat")])
dimnames(geo_dist) <- list(metadata$sample_id, metadata$sample_id)

# test:
for(i in 1:length(metadata$sample_id)){
	if(geo_dist[metadata$sample_id[i],metadata$sample_id[i]] != 0){
		warning("distance between a point and itself is not zero, this indicates a problem")
	}
}


for(i in 1:length(metadata$sample_id)){
	
}


comm_dist <- as.matrix(vegdist(otu_table, method = "bray", diag = TRUE, upper = TRUE))
if(!(identical(dimnames(comm_dist), dimnames(geo_dist)))){
	warning("Whoa there! the row/column names of the two distance matrices do not seem to add up. This is bad.")
}

pdf(file = file.path(fig_dir, "comm_diss_by_geo_dist.pdf"), width = 8, height = 3)
	par(mar = c(4,5,1,1))
	plot(
		x = log(geo_dist + 100), 
		y = comm_dist, 
		ylim = c(0,1), 
		xaxt = "n", 
		las = 1, 
		cex = 0.5, 
		col = rgb(0,0,0,alpha = 0.1), 
		xlab = "Distance between samples (meters)", 
		ylab = "Bray-Curtis dissimilarity"
	)
	axis(side = 1, at = unique(log(metadata$dist_from_shore + 100)), labels = unique(metadata$dist_from_shore))
dev.off()
