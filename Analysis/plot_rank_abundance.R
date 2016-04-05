
my_otu_table <- otu_scaled

threshold_1 <- 100
threshold_2 <- 1000
threshold_3 <- 10000

gt_t1 <- sum(colSums(my_otu_table) > threshold_1)
gt_t2 <- sum(colSums(my_otu_table) > threshold_2)
gt_t3 <- sum(colSums(my_otu_table) > threshold_3)

plot_scale <- 1e6

plot_values <- sort(colSums(my_otu_table)/plot_scale, decreasing = TRUE)

output_path <- file.path(fig_dir, "rank_abundance_rescaled.pdf")
pdf(file = output_path, width = 8, height = 4)
par(mar = c(5, 5, 1, 1))
barplot(
  plot_values, 
  ylim = c(0,3), 
  xaxt = "n", 
  # yaxt = "n", 
  xlab = paste(length(plot_values), "OTUs sorted by decreasing abundance"), 
  ylab = "Millions of reads across samples", 
  las = 1)

# axis(side = 2, at = pretty(range(plot_values)), las = 1) #, lwd = 0, lwd.ticks = 1

text1 <- paste(gt_t1, "sequences >", threshold_1)
text2 <- paste(gt_t2, "sequences >", threshold_2)
text3 <- paste(gt_t3, "sequences >", threshold_3)

legend("topright", legend = c(text1, text2, text3), bty = "n")

text(x = length(plot_values), y = 10000/plot_scale, labels = threshold_3, col = "red", pos = 3)
abline(h = threshold_3/plot_scale, col = "red")

dev.off()