#!/usr/bin/env Rscript

# plot map of sampled site

tif_file <- file.path(data_dir, "ngdc_pug_snd_dm_subset.tif")

library(rgdal) # readGDAL()

tif_import <- readGDAL(tif_file)

class(tif_import)

summary(tif_import)

ignore_above_0 <- function(x){
# convert any values of a vector that are > 0 to 0
  x[ x > 0 ] <- 0
  return(x)
}

attributes(tif_import)$data$band1 <- ignore_above_0(attributes(tif_import)$data$band1)
attributes(tif_import)$data$band1 <- attributes(tif_import)$data$band1/10


library(raster)
data_raster <- raster(tif_import)

# create bounding box
range_lon <- c(-122.48,-122.36)
range_lat <- c(47.68, 47.76)

# create clipping polygon
CP <- extent(c(range_lon, range_lat))

data_raster <- crop(data_raster, CP)

land <- rasterToPolygons(data_raster, fun = function(x) {x >=0 }, dissolve = TRUE)
# shore <- rasterToPolygons(data_raster, fun = function(x) {x ==0 }, dissolve = TRUE)
# class(land)

breakpoints <- seq(from = -325, to = 0, by = 25)

n_col <- length(breakpoints) - 1
library(colorspace)
# col_depth <- sequential_hcl(length(breakpoints) - 1)
pal <- choose_palette()
col_depth <- pal(n_col)
dput(col_depth)
# 13: 
col_depth <- c("#3353B2", "#3E59B3", "#485FB3", "#5165B5", "#5A6CB7", "#6272B8", 
"#6B79BB", "#7380BD", "#7C88C0", "#868FC2", "#9098C5", "#9BA2C9", 
"#AAAFCE")
# 12: col_depth <- c("#294698", "#3A4F9B", "#47599F", "#5463A4", "#606DA9", "#6C77AE", "#7781B4", "#838CBA", "#8F97C0", "#9CA2C6", "#A8AECD", "#B6BAD4")
col_depth <- sequential_hcl(n_col, 
  h = 260, 
  c. = c(60, 30), 
  l = c(20, 70), 
  power = 1, gamma = NULL, fixup = TRUE, alpha = 1)
r <- raster(nrows=1, ncols= n_col)
r <- setValues(r, 1:ncell(r))
plot(r, col = col_depth)

out_file <- file.path(fig_dir, "site_map.pdf")
pdf(file = out_file)
par(mar = c(5,6,1,4))

# this sucks; look into rastervis
plot(x = data_raster_sub, 
  breaks = breakpoints, col = col_depth, 
  xlim = c(-122.46, -122.36),  
  las = 1, xlab = "Latitude", ylab = "Longitude", ann = FALSE,
  axes = TRUE, box = FALSE)

points(metadata[,colname_lon], metadata[,colname_lat], col = "orangered", pch = 4)
plot(land, col = "cornsilk3", add = TRUE)
# plot(shore, col = "black", add = TRUE)

title(xlab = "Latitude")
title(ylab = "Longitude", line = 4)
dev.off()

# plot raster arguments:
# function (x, col, add = FALSE, legend = TRUE, horizontal = FALSE, 
#     legend.shrink = 0.5, legend.width = 0.6, legend.mar = ifelse(horizontal, 
#         3.1, 5.1), legend.lab = NULL, graphics.reset = FALSE, 
#     bigplot = NULL, smallplot = NULL, legend.only = FALSE, lab.breaks = NULL, 
#     axis.args = NULL, legend.args = NULL, interpolate = FALSE, 
#     box = TRUE, breaks = NULL, zlim = NULL, zlimcol = NULL, fun = NULL, 
#     asp, colNA = NA, ...)
