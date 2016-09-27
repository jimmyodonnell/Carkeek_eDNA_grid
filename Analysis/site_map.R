#!/usr/bin/env Rscript

# plot map of sampled site

tif_file <- file.path(data_dir, "ngdc_pug_snd_dm_subset.tif")

library(raster)
data_raster <- raster(tif_file)

ignore_above_0 <- function(x){
# convert any values of a vector that are > 0 to 0
  x[ x > 0 ] <- 0
  return(x)
}

# ignore elevation above sea level
data_raster <- ignore_above_0(data_raster)

# convert from decimeters to meters
data_raster <- data_raster/10

# create bounding box
range_lon <- c(-122.48,-122.36)
range_lat <- c(47.68, 47.76)

# create clipping polygon
CP <- extent(c(range_lon, range_lat))

# crop raster
data_raster <- crop(data_raster, CP)

n_col <- 10
library(colorspace)

col_depth <- sequential_hcl(n_col, 
  h = 260, 
  c. = c(60, 30), 
  l = c(20, 70), 
  power = 1, gamma = NULL, fixup = TRUE, alpha = 1)
col_depth[length(col_depth)] <- "cornsilk3" # set the color of the land

# check the colors
# r <- raster(nrows=1, ncols= n_col)
# r <- setValues(r, 1:ncell(r))
# plot(r, col = col_depth)

# set points to add
mypoints <- SpatialPoints(cbind(metadata[["mean"]][,colname_lon], metadata[["mean"]][,colname_lat]))

#-------------------------------------------------------------------------------
# PLOTTING
plot_name   <- "site_map"

if(!exists("legend_text")){
  legend_text <- list()
}

legend_text[plot_name] <- {
"Map of study area. Depth in meters below sea level is indicated by shading and 25 meter contours. Sampled locations are indicated by red points."
}

EXPORT <- FALSE
if(EXPORT){
  pdf_file    <- file.path(fig_dir, paste(plot_name, ".pdf", sep = ""))
  legend_file <- file.path(fig_dir, paste(plot_name, "_legend.txt", sep = ""))
  writeLines(legend_text[[plot_name]], con = legend_file)
  pdf(file = file.path(fig_dir, plot_pdf), width = 6, height = 5) #, width = 8, height = 3  
}
par(mar = c(5,6,1,4))

library(rasterVis)
myTheme <- rasterTheme(region = col_depth)
# myTheme$add.line$col <- "pink" # contour line colors
levelplot(
  x = data_raster, 
  margin = FALSE, contour = TRUE, 
  par.settings = myTheme
) + 
layer(sp.points(mypoints, col = "orangered", pch = 4))

if(EXPORT){
  dev.off()
}

# plot raster arguments:
# function (x, col, add = FALSE, legend = TRUE, horizontal = FALSE, 
#     legend.shrink = 0.5, legend.width = 0.6, legend.mar = ifelse(horizontal, 
#         3.1, 5.1), legend.lab = NULL, graphics.reset = FALSE, 
#     bigplot = NULL, smallplot = NULL, legend.only = FALSE, lab.breaks = NULL, 
#     axis.args = NULL, legend.args = NULL, interpolate = FALSE, 
#     box = TRUE, breaks = NULL, zlim = NULL, zlimcol = NULL, fun = NULL, 
#     asp, colNA = NA, ...)
