# calculate Moran's I

my_metadata <- metadata[[dataset]]

# grab diversity data
temp <- do.call("cbind.data.frame",  lapply(div_metrics, function(x) cbind.data.frame(x[5:length(x)])))

# bind it all together
my_data <- data.frame(
  site = my_metadata[, colname_env_sample], 
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

# calculate Moran's I
library(ape)
moran_out <- apply(
  my_data[,4:ncol(my_data)], 2, 
  function(x) Moran.I(x, dist_sphere_inv)
)

# to see all p-values
sapply(moran_out, "[[", "p.value")

