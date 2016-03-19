# INPUT COMES FROM SCRIPT "load_data.r"

# plot scaled proportional abundance in space (i.e. darkest color is wherever OTU's proportional abundance is at maximum)

for(i in c(1:5,7:10)){
	plot(
		metadata$lon, 
		metadata$lat, 
		bg = rgb(r = 1, g = 0, b = 0, alpha = otu_table_prop[metadata$sample_id, i]/max(otu_table_prop[metadata$sample_id, i])), 
		pch = 21, 
		main = colnames(otu_table_prop)[i], 
		xlab = "Longitude", 
		ylab = "Latitude"
	)
}

otu_table[,6]

pdf(file = file.path(fig_dir, "otu_in_space.pdf"), width = 7, height = 7)
	plot(
		metadata$lon, 
		metadata$lat, 
		bg = rgb(r = 1, g = 0, b = 0, alpha = otu_table_prop[metadata$sample_id, "DUP_3"]/max(otu_table_prop[metadata$sample_id, "DUP_3"])), 
		pch = 21, 
		cex = otu_table_prop[metadata$sample_id, "DUP_3"]/max(otu_table_prop[metadata$sample_id, "DUP_3"])*1.5, 
		# cex = 1.5, 
		main = "DUP_3", 
		xlab = "Longitude", 
		ylab = "Latitude"
	)
	points(
		metadata$lon, 
		metadata$lat, 		
		pch = ifelse(otu_table_prop[metadata$sample_id, "DUP_3"]/max(otu_table_prop[metadata$sample_id, "DUP_3"]) == 0, 4, NA_integer_) 
	)
dev.off()

set.seed(8) # to control jitter
pdf(file = file.path(fig_dir, "otu_in_space.pdf"), width = 7, height = 7)
	plot(
		jitter(x = metadata$dist_along_shore, factor = 0.2), 
		xaxt = "n", 
		xlim = c(-500, 2500), 
		y = jitter(log(metadata$dist_from_shore + 10), factor = 0.2), 
		yaxt = "n",
		# log = "y",
		bg = rgb(r = 1, g = 0, b = 0, alpha = otu_table_prop[metadata$sample_id, "DUP_3"]/max(otu_table_prop[metadata$sample_id, "DUP_3"])), 
		pch = 21, 
		cex = otu_table_prop[metadata$sample_id, "DUP_3"]/max(otu_table_prop[metadata$sample_id, "DUP_3"])*1.5, 
		# cex = 1.5, 
		main = "DUP_3", 
		las = 1, 
		xlab = "Position along shore (meters)", 
		ylab = "Position from 0 (meters)"
	)
	points(
		x = jitter(metadata$dist_along_shore, factor = 0.2), 
		y = jitter(log(metadata$dist_from_shore + 10), factor = 0.2), 
		# log = "y",
		cex = 0.5, 
		pch = ifelse(otu_table_prop[metadata$sample_id, "DUP_3"]/max(otu_table_prop[metadata$sample_id, "DUP_3"]) == 0, 4, NA_integer_) 
	)
	axis(side = 2, at = unique(log(metadata$dist_from_shore + 10)), labels = unique(metadata$dist_from_shore), las = 1)
	axis(side = 1, at = unique(metadata$dist_along_shore), labels = unique(metadata$dist_along_shore), las = 1)
dev.off()

