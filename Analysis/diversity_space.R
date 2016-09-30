#-------------------------------------------------------------------------------
# Diversity in Space
#-------------------------------------------------------------------------------
EXPORT <- FALSE

dataset <- "mean"

my_metadata <- metadata[[dataset]]

# grab diversity data
temp <- do.call(
  "cbind.data.frame",  
  lapply(
    div_metrics, 
    function(x) {
      cbind.data.frame(x[5:length(x)])
    } # select only datasets of same length
  )
)

# bind it all together
my_data <- data.frame(
  site = my_metadata[, colname_env_sample], 
  dist = my_metadata[, colname_ycoord], 
  lat  = my_metadata[, colname_lat], 
  lon  = my_metadata[, colname_lon]  
)
my_data <- cbind.data.frame(
  my_data, 
  temp[match(rownames(temp), my_data$site),], 
  row.names = NULL
)
rm(temp)

# compute distances
library(geosphere)
dist_sphere <- distm(
  x = my_data[,c("lon", "lat")],
  fun = distHaversine
)
dimnames(dist_sphere) <- list(my_data$site, my_data$site)

# calculate the inverse distances and set diagonal to 0
dist_sphere_inv <- 1/dist_sphere
diag(dist_sphere_inv) <- 0

#-------------------------------------------------------------------------------
# calculate Moran's I
library(ape)
moran_out <- apply(
  my_data[,5:ncol(my_data)], 2, 
  function(x) Moran.I(x, dist_sphere_inv)
)

# to see all p-values
sapply(moran_out, "[[", "p.value")
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# model and plot diversity as a function of distance from shore (full)

div_dist <- list()

plot_name <- "diversity_distance_all"

if(!exists("legend_text")){ legend_text <- list()}
legend_text[plot_name] <- {
"Legend text goes here."
}
if(EXPORT){
  pdf_file    <- file.path(fig_dir, paste(plot_name, ".pdf", sep = ""))
  legend_file <- file.path(fig_dir, paste(plot_name, "_legend.txt", sep = ""))
  writeLines(legend_text[[plot_name]], con = legend_file)
  pdf(file = pdf_file, width = 16, height = 8)
}
par(mfrow = c(3, 6))
for(i in 5:22){
par(mar = c(4,4,1,1))
div_dist[[colnames(my_data)[i]]] <- lm(
  my_data[,i] ~ my_data$dist
)
plot(x = my_data$dist, y = my_data[,i],
  # log = "x", 
  xlab = "Distance from shore (meters)", 
  ylab = colnames(my_data)[i]
)
abline(div_dist[[colnames(my_data)[i]]], col = "red")
}
if(EXPORT){
  dev.off()
}
dataset <- grep(".mean$", names(div_dist))
(summaries <- lapply(div_dist[dataset], summary))

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# model and plot diversity as a function of distance from shore, but
# focus on only a single dataset ("mean")

col_sub <- grep(".mean$", colnames(my_data))

div_dist <- list()

plot_name <- "diversity_distance"

if(!exists("legend_text")){ legend_text <- list()}
legend_text[plot_name] <- {
"Legend text goes here."
}
if(EXPORT){
  pdf_file    <- file.path(fig_dir, paste(plot_name, ".pdf", sep = ""))
  legend_file <- file.path(fig_dir, paste(plot_name, "_legend.txt", sep = ""))
  writeLines(legend_text[[plot_name]], con = legend_file)
  pdf(file = pdf_file, width = 8, height = 2.5)
}

par(mfrow = c(1, length(col_sub)))

for(i in 1:length(col_sub)){
par(mar = c(4,4,1,1))
temp_dat <- data.frame(x = my_data$dist, y = my_data[,col_sub[i]])
metric_name <- strsplit(colnames(my_data)[col_sub[i]], "\\.")[[1]][1]

div_dist[[metric_name]] <- list()
div_dist[[metric_name]]$out <- lm(
  y ~ log(x + 1), data = temp_dat
)
plot(x = temp_dat$x, y = temp_dat$y,
  log = "x", 
  xlab = "Distance from shore (meters)", 
  ylab = metric_name
)

x_pred <- seq(from = 0, to = 5000, length = 101)

div_dist[[metric_name]]$conf <- predict(
    object   = div_dist[[metric_name]]$out, 
    newdata  = data.frame(x = x_pred), 
    interval = "confidence", 
    level    = 0.95
)

plot_model(div_dist[[metric_name]], x_pred, line_type = 2)
# abline(div_dist[[colnames(my_data[,col_sub])[i]]], col = "red")
}
if(EXPORT){
  dev.off()
}
