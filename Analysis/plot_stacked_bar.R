# plot stacked bar

library(ggplot2)
library(reshape2)


mytable <- otu_filt[,1:20]
mytable_long <- melt(mytable)
colnames(mytable_long) <- c("sample_id", "OTU", "count")

plot_dat <- mytable_long


outfile_name <- "stacked_bar_community.pdf"
# pdf(file = file.path(fig_dir, outfile_name), width = 9, height = 5)
p <- ggplot(aes(x = sample_id, weight = count, fill = OTU), data = plot_dat)
p + geom_bar() + 
  scale_fill_discrete("Taxa") + 
  labs(x="Water Samples", y="Number of Sequences") +
  # Remove axis ticks and tick mark labels
  theme(
      axis.text.x = element_blank(),
      # axis.text.y = element_blank(),
      # axis.ticks = element_blank(), 
      panel.background = element_rect(fill = 'white', colour = 'white')
      )
# dev.off()

# note high value for sample "lib_B_tag_GCGCTC" in PCT-C-0500


# head(data_full_long)
# tmp <- data_full_long[order(data_full_long$env_sample_name),]
# head(plot_dat)

# plot_dat <- tmp[tmp$env_sample_name == "PCT-C-0500",]
# plot_dat <- plot_dat[plot_dat$count > 200,]
