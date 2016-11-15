library(vegan)

# requires a list of OTU tables and a list of corresponding metadata

EXPORT <- TRUE

div_metrics <- list()

# div_metrics[["Shannon"]] <- lapply(otu_table, diversity, index = "shannon") 
div_metrics[["Simpson"]] <- lapply(otu_table, function(x){
  1 - diversity(x, index = "simpson")})
div_metrics[["Richness"]] <- lapply(otu_table, function(x) {rowSums(x > 0)})

which_data <- "mean"

plot_name <- "diversity"
if(!exists("legend_text")){legend_text <- list()}
legend_text[plot_name] <- {"Aggregate measures of diversity at each sample site."}

if(EXPORT){
  pdf_file    <- file.path(fig_dir, paste(plot_name, ".pdf", sep = ""))
  legend_file <- file.path(fig_dir, paste(plot_name, "_legend.txt", sep = ""))
  writeLines(legend_text[[plot_name]], con = legend_file)
  pdf(file = pdf_file, width = 5, height = 4) #, width = 8, height = 3
}

the_layout <- layout(mat = matrix(c(1,2), nrow = 1), widths = c(11,6))

margins_subplot <- list(c(4,5,3,5), c(4,0,3,3), c(4,0,3,3))

for(i in 1:length(div_metrics)){

plot_dat <- div_metrics[[i]][[which_data]][
  metadata[[which_data]][,colname_env_sample]
]

n_colors <- 10
#Create a function to generate a continuous color palette
if(require(viridis)){
  colorpalette <- viridis
} else {
  colorpalette <- colorRampPalette(c('darkorchid','aquamarine'))
}

colorified <- colorpalette(n_colors)[cut(plot_dat, breaks = n_colors)]

y_for_plotting <- as.numeric(as.factor(metadata[[which_data]][ , colname_ycoord]))

par(mar = margins_subplot[[i]], xpd = TRUE)
plot(
  x = metadata[[which_data]][,colname_xcoord], 
  xlim = c(-500, 2500), # add a little space to the right and left
  y = y_for_plotting, 
  col = colorified,
  pch = 19, 
  cex = 3, 
  xlab = "", 
  ylab = "", 
  axes = FALSE, 
  las = 1
)
if(i == 1){
  axis(side = 2, at = unique(y_for_plotting), 
    labels = unique(metadata[[which_data]][,colname_ycoord]), las = 1)
  title(ylab = "Distance from 0 (meters)")
}

title(main = names(div_metrics[i]))
title(xlab = "Transect", line = 1)

# add color legend
dat_col_leg <- seq(min(plot_dat), max(plot_dat), length.out = n_colors)
colvec <- colorpalette(n_colors)[cut(dat_col_leg, breaks = n_colors)]
color_legend(col = colvec, lev = sort(plot_dat))


}
if(EXPORT){
  dev.off()
}
