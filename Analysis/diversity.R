library(vegan)

# requires a list of OTU tables and a list of corresponding metadata

shannon <- lapply(otu_table, diversity, index = "shannon") 
simpson <- lapply(otu_table, diversity, index = "simpson") 
# "shannon", "simpson", "invsimpson"

richness <- lapply(otu_table, function(x) {rowSums(x > 0)})


which_data <- "mean_unfilt"

plot_dat <- simpson[[which_data]][
  metadata[[which_data]][,colname_env_sample]
]

n_colors <- 10
#Create a function to generate a continuous color palette
colorpalette <- colorRampPalette(c('blue','red'))

colorified <- colorpalette(n_colors)[cut(plot_dat, breaks = n_colors)]

par(mar = c(4,5,1,3))
plot(
  x = 
metadata[[which_data]][,colname_xcoord], 
  y = 
metadata[[which_data]][,colname_ycoord], 
  col = colorified,
  pch = 19, 
  cex = 3, 
  xlab = "Transect", 
  ylab = "Distance from 0", 
  axes = FALSE, 
  las = 1
)
axis(2, las = 1)

# add color legend
colvec <- colorpalette(n_colors)[cut(sort(plot_dat), breaks = n_colors)]
color_legend(col = colvec, lev = sort(plot_dat))