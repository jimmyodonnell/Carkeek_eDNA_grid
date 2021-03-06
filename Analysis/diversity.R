library(vegan)

# requires a list of OTU tables and a list of corresponding metadata

EXPORT <- FALSE

div_metrics <- list()

# div_metrics[["Shannon"]] <- lapply(otu_table, diversity, index = "shannon") 
div_metrics[["Evenness"]] <- lapply(otu_table, function(x){
  1 - diversity(x, index = "simpson")})
div_metrics[["Richness"]] <- lapply(otu_table, function(x) {rowSums(x > 0)})

which_data <- "mean"

plot_name <- "diversity"
if(!exists("legend_text")){legend_text <- list()}
legend_text[plot_name] <- {
"Aggregate measures of diversity at each sample site.
Data are rarefied counts of mitochondrial 16S sequences collected from 3 parallel transects in Puget Sound, Washington, USA.
Evenness (A) is the probability that two sequences drawn at random are different; richness (B) represents the total number of unique sequences from that location."
}

if(EXPORT){
  pdf_file    <- file.path(fig_dir, paste(plot_name, ".pdf", sep = ""))
  legend_file <- file.path(fig_dir, paste(plot_name, "_legend.txt", sep = ""))
  writeLines(legend_text[[plot_name]], con = legend_file)
  pdf(file = pdf_file, width = 5, height = 4) #, width = 8, height = 3
}

the_layout <- layout(mat = matrix(c(1,2), nrow = 1), widths = c(11,6))

margins_subplot <- list(c(3,5,1,5), c(3,0,1,3))

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

# add letter for legend
mtext(
  text = LETTERS[i], 
  side = 1, 
  adj = 1.4, 
  # at = c(max(metadata[[which_data]][,colname_xcoord])*1.2, 0), 
  line = 1, 
  cex = 2 )

title(xlab = "Transect", line = 1)

# add color legend
dat_col_leg <- seq(min(plot_dat), max(plot_dat), length.out = n_colors)
colvec <- colorpalette(n_colors)[cut(dat_col_leg, breaks = n_colors)]
color_legend(col = colvec, lev = sort(plot_dat))


}
if(EXPORT){
  dev.off()
}
