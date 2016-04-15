my_table    <- otu_mean
my_metadata <- metadata_mean

library(geoR) # variog

pdf(file = file.path(fig_dir, "variograms_taxa.pdf"))
for(i in 1:20){

	vg <- variog(coords = my_metadata[,c(colname_lon, colname_lat)], data = my_table[,i])

	plot(vg, type = "b", main = colnames(my_table)[i])

}
dev.off()
# vg.summary <- cbind(c(1:10), v2$v, v2$n)
# colnames(v2.summary) <- c("lag", "semi-variance", "# of pairs")
# v2.summary

# plot(v2, type = "b", main = "Variogram: Av8top")
