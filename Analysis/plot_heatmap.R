# Heat Mapping in ggplot, courtesy of Ole Shelton
#-------------------------------------------------------------------------------

library(ggplot2) # ggplot()
library(reshape2) # melt()

dis_bray <- vegdist(otu_filt, method = "bray", binary = FALSE, diag = TRUE, upper = TRUE)
plot_data <- melt(as.matrix(dis_bray))
head(plot_data)
names(plot_data)[3] <- "dissimilarity"

pdf(file = file.path(fig_dir,"heatmap_sites.pdf"))
p.all <- ggplot(plot_data,
                 aes(x = Var1, y = Var2, fill = dissimilarity)) +
  geom_tile() +
  scale_fill_gradientn(colours = c("blue",grey(0.8),"red"),limits=c(0,1),na.value="black") +
    coord_equal() +
  # scale_x_discrete(breaks=1:11)+
  # scale_y_discrete(breaks=1:11,limits = rev(levels(plot_data$Var1))) +
  labs(x= "Site", y= "Site") +
  ggtitle("Site Dissimilarity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

print(p.all)
dev.off()
