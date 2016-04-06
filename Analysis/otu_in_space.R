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

pdf(file = file.path(fig_dir, "otu_in_space_mean.pdf"), width = 7, height = 7)
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

################################################################################
# LOOP TO PRINT MANY AT THE SAME TIME
my_otus <- otu_mean
rownames(otu_mean) # PCT-C-0000 etc, aka "env_sample_name"
my_metadata <- metadata_exp[!duplicated(metadata_exp[,"env_sample_name"]),]

pdf(file = file.path(fig_dir, "otu_in_space_50.pdf"), width = 7, height = 7)
for(i in 1:50){
	plot(
		x = my_metadata$dist_along_shore, 
		xaxt = "n", 
		xlim = c(-500, 2500), 
		y = log(my_metadata$dist_from_shore + 10), 
		yaxt = "n",
		# log = "y",
		bg = rgb(r = 1, g = 0, b = 0, alpha = my_otus[my_metadata$env_sample_name, i]/max(my_otus[my_metadata$env_sample_name, i])), 
		pch = 21, 
		cex = my_otus[my_metadata$env_sample_name, i]/max(my_otus[my_metadata$env_sample_name, i])*5, 
		# cex = 1.5, 
		main = colnames(my_otus)[i], 
		las = 1, 
		xlab = "Position along shore (meters)", 
		ylab = "Position from 0 (meters)"
	)
	points(
		x = my_metadata$dist_along_shore, 
		y = log(my_metadata$dist_from_shore + 10), 
		# log = "y",
		cex = 1, 
		pch = ifelse(
			my_otus[my_metadata$env_sample_name, i]/
				max(my_otus[my_metadata$env_sample_name, i]) == 0, 
			4, NA_integer_) 
	)
	axis(side = 2, at = unique(log(my_metadata$dist_from_shore + 10)), labels = unique(my_metadata$dist_from_shore), las = 1)
	axis(side = 1, at = unique(my_metadata$dist_along_shore), labels = unique(my_metadata$dist_along_shore), las = 1)
}
dev.off()

