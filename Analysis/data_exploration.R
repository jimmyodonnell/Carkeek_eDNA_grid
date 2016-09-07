
REMOVE_HUMAN_DNA <- FALSE
if(REMOVE_HUMAN_DNA){
	mytable <- as.matrix(otu_table_raw[,-1])# no human
} else {
	mytable <- as.matrix(otu_table_raw) # [,1:200]
}

if(!identical(rownames(mytable), metadata[,colname_sampleid])){
  warning("the rownames of the table and the sample ID column should be in the same order")
}
order_by_env <- unlist(split(rownames(mytable), metadata[,colname_env_sample]))
mytable <- mytable[order_by_env,]
rownames(mytable) <- names(order_by_env)
seq_depth <- rowSums(mytable)
seq_depth_df <- data.frame(sample = names(seq_depth), seq_depth = seq_depth/max(seq_depth), row.names = NULL)

mytable <- mytable/rowSums(mytable)
other <- rowSums(mytable[,20:ncol(mytable)])
mytable <- mytable[,1:19]
mytable <- cbind(mytable, other)

lines_env_samples <- which(!duplicated(metadata[ , colname_env_sample])) - 0.5

## Plot Sequencing depth
pdf(file = file.path(fig_dir, "sequencing_depth.pdf"), width = 5, height = 20)
par(mar = c(5, 8, 1, 1))
stripchart(x = as.list(seq_depth), las = 1, xlab = "number of sequences", axes = FALSE, pch = 1) # , axes = FALSE , cex.axis = 0.6
axis(side = 1)
axis(side = 2, at = seq_along(seq_depth), labels = names(seq_depth), las = 1, cex.axis = 0.7)
abline(v = seq(from = 0, to = max(seq_depth), by = 100000), col = hsv(1,0.1,0.1, alpha = 0.1))
abline(h = lines_env_samples, col = hsv(1,0.1,0.1, alpha = 0.1))
dev.off()

library(reshape2)
mytable_long <- melt(mytable)
colnames(mytable_long) <- c("sample_id", "OTU", "count")

plot_dat <- mytable_long

gghue <- function(n){
	hues = seq(15, 375, length = n+1)
	hcl(h = hues, l = 65, c = 100)[1:n]
}

otu_cols <- gghue(ncol(mytable)-1)
otu_cols <- c(otu_cols, "gray")
library(ggplot2)
if(REMOVE_HUMAN_DNA){
  outfile_name <- "stacked_bar_community_horiz_nohuman.pdf"
} else {
  outfile_name <- "stacked_bar_community_horiz.pdf"
}
pdf(file = file.path(fig_dir, outfile_name), width = 5, height = 20)
p <- ggplot() +
  geom_bar(data = plot_dat, aes(x = sample_id, weight = count, fill = OTU)) +
# p + geom_bar() + 
  scale_fill_manual(values = otu_cols) + 
  # scale_fill_discrete("Taxa") + 
  coord_flip() +
  labs(x="Bars = PCR replicates; lines demark DNA extracts; x = sequencing depth", y="Proportion of Sequences") +
  # Remove axis ticks and tick mark labels
  theme(
      # axis.text.x = element_blank(),
      # axis.text.y = element_blank(),
      # axis.ticks = element_blank(), 
      panel.background = element_rect(fill = 'white', colour = 'white')
      ) + 
  geom_vline(xintercept = lines_env_samples, colour = hsv(0,0,0.5)) + 
  geom_point(data = seq_depth_df, aes(x = sample, y = seq_depth), shape = 4)
p
dev.off()

