#!/usr/bin/env Rscript
#-------------------------------------------------------------------------------

the_otu_table <- otu_filt
the_metadata  <- metadata_mean


#-------------------------------------------------------------------------------
data_by_dist <- split(
	x = as.data.frame(the_otu_table), 
	f = metadata_mean[,colname_ycoord]
)

distances <- rep(names(data_by_dist), times = sapply(data_by_dist, nrow))

dis_by_dist <- lapply(data_by_dist, function(x) 
	as.vector(vegdist(x, method = "bray"))
)

to_plot <- data.frame(dist = as.numeric(distances), vegdist = unlist(dis_by_dist))

# transform distance to shore 
lm_out <- lm(to_plot$vegdist ~ log(to_plot$dist + 1))
summary(lm_out)
# no significant slope


pdf(file = file.path(fig_dir, "dissimilarity_from_shore.pdf"),
	width = 4, 
	height = 4
)
par(mar = c(4,4,1,1))
plot(
	x = to_plot$dist,
	xlab = "Distance from shore", 
	y = to_plot$vegdist, 
	ylab = "Bray-Curtis Dissimilarity", 
	ylim = c(0,1), 
	pch = 21, 
	col = hsv(1,1,0), 
	bg  = hsv(1,1,0,0.2), 
	# log = "x"
	las = 1
)
points(
	y = sapply(dis_by_dist, mean), 
	x = as.numeric(names(dis_by_dist)), 
	pch = 23, 
	col = hsv(1,1,1), 
	bg  = hsv(1,1,1,0.2)
)
legend(
	"bottomright", 
	bty = "n", 
	legend = "mean", 
	pch = 23, 
	col = hsv(1,1,1), 
	pt.bg  = hsv(1,1,1,0.2)
)
# abline(lm_out, col = hsv(0.6, 1, 1), lwd = 2, lty = 2)
dev.off()

#-------------------------------------------------------------------------------
# individual otu abundance variance
var_by_dist <- lapply(data_by_dist, function(x)
	apply(x, MARGIN = 2, FUN = function(x) var(log(x + 1)))
)

distances <- rep(names(var_by_dist), times = sapply(var_by_dist, length))

to_plot <- data.frame(dist = as.numeric(distances), var = unlist(var_by_dist))

plot(
	x = to_plot$dist,
	xlab = "Distance from shore", 
	y = to_plot$var, 
	ylab = "Within-OTU Variance of log(counts) (mean = red diamond)", 
	# ylim = c(0,1), 
	# log = "x"
	pch = 21, 
	col = hsv(1,1,0), 
	bg  = hsv(1,1,0,0.2), 
	las = 1
)
points(
	y = sapply(var_by_dist, mean), 
	x = as.numeric(names(var_by_dist)), 
	pch = 23, 
	col = hsv(1,1,1), 
	bg  = hsv(1,1,1,0.2)
)

# abline(lm_out, col = hsv(0.6, 1, 1), lwd = 2, lty = 2)

pdf(file = file.path(fig_dir, "variance_from_shore.pdf"),
	width = 4, 
	height = 4
)
par(mar = c(4,4,1,1))
stripchart(
	x = var_by_dist, 
	ylab = "Within-OTU Variance of counts log(x+1)", 
	xlab = "Distance from shore", 
	method = "jitter", 
	pch = 21, 
	col = hsv(1,1,0), 
	bg  = hsv(1,1,0,0.2), 
	las = 1, 
	vertical = TRUE
)
points(
	x = sapply(var_by_dist, mean), 
	pch = 23, 
	col = hsv(1,1,1), 
	bg  = hsv(1,1,1,0.2)
)
legend(
	"topleft", 
	bty = "n", 
	legend = "mean", 
	pch = 23, 
	col = hsv(1,1,1), 
	pt.bg  = hsv(1,1,1,0.2)
)
dev.off()
