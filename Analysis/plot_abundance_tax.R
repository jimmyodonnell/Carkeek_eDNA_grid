# plot sorted abundance

plot_dat <- colSums(otu_mean_LCT[,1:20])
outfile_name <- "rank_abundance_tax20.pdf"
pdf(file = file.path(fig_dir, outfile_name))
par(mar = c(4, 10, 1, 1))
barplot(
	plot_dat,
	horiz = TRUE,
	names.arg = strtrim(names(plot_dat), width = 18),
	xlab = "number of sequences",
	las = 1
	)
dev.off()
