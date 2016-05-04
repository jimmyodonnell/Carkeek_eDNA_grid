my_table    <- otu_filt
my_metadata <- metadata_mean

library(geoR) # variog

# Attempting multivariate variogram
# TODO as of now, I think it just calculates individual variograms
my_geodata <- as.geodata(
	obj = cbind.data.frame(
		my_metadata[,c(colname_lon, colname_lat)], #c(colname_xcoord, colname_ycoord)
		my_table), 
	coords.col = 1:2, 
	data.col = 2 + 1:ncol(my_table)
)
vg_all <- variog(my_geodata)
plot(vg_all, type = "b")

# variogram for individual taxon abundance
pdf(file = file.path(fig_dir, "variogram_taxa.pdf"))
for(i in 1:20){

	vg <- variog(coords = my_metadata[,c(colname_lon, colname_lat)], data = my_table[,i])

	plot(vg, type = "b", main = colnames(my_table)[i])

}
dev.off()


# variogram for shannon index
pdf(file = file.path(fig_dir, "variogram_shannon.pdf"))
shannon <- diversity(my_table, index = "shannon")
vg <- variog(coords = my_metadata[,c(colname_lon, colname_lat)], data = my_table[,i])
plot(vg, type = "b", main = "Variogram of Shannon Index")
dev.off()


# variogram for taxon richness (number of taxa w/ mean abundance > 1)
thresholds <- c(0, 1, 10, 100, 1000, 10000)
pdf(file = file.path(fig_dir, "variogram_richness.pdf"))
for(i in 1:length(thresholds)){
	richness <- rowSums(my_table > i)
	vg <- variog(coords = my_metadata[,c(colname_lon, colname_lat)], data = richness)
	plot(vg, type = "b", main = paste("Variogram of Richness, abundance threshold = ", thresholds[i]))
}
dev.off()
