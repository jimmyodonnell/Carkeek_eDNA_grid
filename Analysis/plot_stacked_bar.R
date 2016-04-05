# plot stacked bar

head(data_full_long)


tmp <- data_full_long[order(data_full_long$env_sample_name),]
head(plot_dat)
plot_dat <- tmp[tmp$env_sample_name == "PCT-C-0500",]
plot_dat <- plot_dat[plot_dat$count > 200,]
library(ggplot2)
p <- ggplot(aes(x = sample_id, weight = count, fill = OTU), data = plot_dat)
p + geom_bar() + 
  scale_fill_discrete("Legend Title") + 
  labs(x="X Label", y="Y Label", title="An Example Stacked Column Chart")

# note high value for sample "lib_B_tag_GCGCTC" in PCT-C-0500
