# dissimilarity by distance for reduced data set (mean abundance within sample for each OTU)

library(geosphere) # distm()

my_otu <- otu_mean
rownames(my_otu) # should be e.g. PCT-C-0000 etc, aka "env_sample_name"
my_metadata <- metadata_mean# metadata[!duplicated(metadata[,colname_env_sample]),]

# calculate pairwise great circle distance between sampling locations using Haversine method
geo_dist <- as.dist(distm(x = my_metadata[,c(colname_lon, colname_lat)], fun = distHaversine))
# dimnames(geo_dist) <- list(my_metadata$env_sample_name, my_metadata$env_sample_name)

comm_dist <- vegdist(my_otu, method = "bray") #, diag = TRUE, upper = TRUE
if(!(identical(dimnames(comm_dist), dimnames(geo_dist)))){
	warning("Whoa there! the row/column names of the two distance matrices do not seem to add up. This is bad.")
}

# fit michaelis menten curve
the_data <- data.frame(comm = as.vector(comm_dist), space = as.vector(geo_dist))
mm_fit <- nls(
  formula = comm ~ Vm * space/(Km + space),
  data = the_data,
  start = list(
    Vm = max(the_data$comm), # comment out and set Vm to 1 if asymptote should be 1
    Km = max(the_data$comm)/2
  )
)
mm_prediction <- predict(mm_fit)

geo_dist_scaled <- log(geo_dist + 100)
plot_x <- geo_dist # geo_dist_scaled

pdf(file = file.path(fig_dir, "diss_by_dist.pdf")) #, width = 8, height = 3
	par(mar = c(4,5,1,1))
	plot(
		x = plot_x,
		y = comm_dist,
		ylim = c(0,1),
		xaxt = "n",
		pch = 21,
		las = 1,
		cex = 1,
		col = 1,
		bg = rgb(0,0,0,alpha = 0.1 ), #,alpha = 0.1
		xlab = "Distance between samples (meters)",
		ylab = "Bray-Curtis dissimilarity"
	)
	axis(side = 1)
	#, at = unique(log(metadata$dist_from_shore + 100)), labels = unique(metadata$dist_from_shore)
# abline(v = unique(log(metadata$dist_from_shore + 100)))

# Add Michaelis Menten Fit
lines(x = c(0, sort(geo_dist)), y = c(0, sort(mm_prediction)), col = "red", lwd = 2, lty = 2)

# Add LOESS line
# smoothed <- loess.smooth(
				# x = plot_x,
				# y = comm_dist,
				# span = 2/3,
				# degree = 1,
				# family = "gaussian"
				# )
# lines(smoothed, col="red", lwd=2)
dev.off()

#-------------------------------------------------------------------------------
# Some Possible analyses
#-------------------------------------------------------------------------------
lm_out <- lm(log(1/comm_dist) ~ geo_dist)
# Then the regression coefficient is usually used in the literature as the descriptor of distance decay, or the distance at which 50% of the maximum similarity is observed.
summary(lm_out)

mantel(comm_dist, geo_dist, perm = 9999)


#-------------------------------------------------------------------------------
# SPLIT BY TRANSECT
#-------------------------------------------------------------------------------
transects <- c("PCT-S", "PCT-C", "PCT-N")

diss_by_transect <- list()
for(transect in 1:length(transects)){
	diss_by_transect[[transect]] <- vegdist(my_otu[grep(transects[transect], rownames(comm_dist)),], method = "bray")
}

dist_by_transect <- list()
for(transect in 1:length(transects)){
	transect_coords <- my_metadata[
		grep(transects[transect], my_metadata[, colname_env_sample]),
		c(colname_lon, colname_lat)]
	dist_by_transect[[transect]] <- as.dist(distm(x = transect_coords, fun = distHaversine))
}

# function for plot colors (sorta like ggplot)
gghue <- function(n){
	hues = seq(15, 375, length = n+1)
	hcl(h = hues, l = 65, c = 100)[1:n]
}
mycols <- gghue(length(transects))

# mycols <- c("rgb(1,0,0)", "rgb(0,1,0)", "rgb(0,0,1)")

# pdf(file = file.path(fig_dir, "diss_by_dist_by_transect.pdf"), width = 8, height = 3)
par(mar = c(4,5,1,1))
	plot(
		x = geo_dist,
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
	axis(side = 1) #, at = unique(log(metadata$dist_from_shore + 100)), labels = unique(metadata$dist_from_shore)

for( i in 1:length(transects)){
	points(
		x = dist_by_transect[[i]],
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
legend("bottomright", legend = transects, bty = "n", pch = 19, col = mycols)

smoothed <- list()
for(i in 1:length(transects)) {
	smoothed[[i]] <- loess.smooth(
				x = dist_by_transect[[i]],
				y = diss_by_transect[[i]],
				span = 1,
				degree = 1,
				family = "gaussian"
				)
	lines(smoothed[[i]], col = mycols[[i]], lwd = 2)
}
# dev.off()
