#-------------------------------------------------------------------------------
# SPLIT BY TRANSECT
#-------------------------------------------------------------------------------
split_by_transect <- FALSE
if(split_by_transect){
  transects <- c("PCT-S", "PCT-C", "PCT-N")
  
  diss_by_transect <- list()
  for(transect in 1:length(transects)){
    diss_by_transect[[transect]] <- vegdist(
      my_table[grep(transects[transect], rownames(as.matrix(comm_dist))),], 
      method = "bray")
  }
  
  dist_by_transect <- list()
  for(transect in 1:length(transects)){
    transect_coords <- my_metadata[
      grep(transects[transect], my_metadata[, colname_env_sample]),
      c(colname_lon, colname_lat)]
    dist_by_transect[[transect]] <- as.dist(distm(x = transect_coords, fun = distHaversine))
  }
  
  # function for plot colors (sorta like ggplot)
  gghue <- function(n){
    hues = seq(15, 375, length = n+1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  mycols <- gghue(length(transects))
  
  # mycols <- c("rgb(1,0,0)", "rgb(0,1,0)", "rgb(0,0,1)")
  
  if(export_plots){
    pdf(file = file.path(fig_dir, "diss_by_dist_by_transect.pdf"), width = 8, height = 3)
  }
  
  par(mar = c(4,5,1,1))
  plot(
    x = geo_dist,
    y = comm_dist,
    ylim = c(0,1),
    xaxt = "n",
    pch = "",
    las = 1,
    cex = 1,
    col = rgb(0,0,0), #,alpha = 0.1
    xlab = "Distance between samples (meters)",
    ylab = "Bray-Curtis dissimilarity"
  )
  axis(side = 1) #, at = unique(log(metadata$dist_from_shore + 100)), labels = unique(metadata$dist_from_shore)
  
  for( i in 1:length(transects)){
    points(
      x = dist_by_transect[[i]],
      y = diss_by_transect[[i]],
      ylim = c(0,1),
      xaxt = "n",
      las = 1,
      pch = 19,
      cex = 1,
      col = mycols[i], #,alpha = 0.1
      xlab = "Distance between samples (meters)",
      ylab = "Bray-Curtis dissimilarity"
    )
  }
  legend("bottomright", legend = transects, bty = "n", pch = 19, col = mycols)
  
  smoothed <- list()
  for(i in 1:length(transects)) {
    smoothed[[i]] <- loess.smooth(
      x = dist_by_transect[[i]],
      y = diss_by_transect[[i]],
      span = 1,
      degree = 1,
      family = "gaussian"
    )
    lines(smoothed[[i]], col = mycols[[i]], lwd = 2)
  }
  if(export_plots){
    dev.off()
  }
  
}
