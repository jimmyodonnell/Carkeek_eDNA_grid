my_table    <- otu_mean
my_metadata <- metadata_mean

library(geoR) # variog

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
