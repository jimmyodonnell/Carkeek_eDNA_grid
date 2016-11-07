#-------------------------------------------------------------------------------
# Diversity in Space
#-------------------------------------------------------------------------------
EXPORT <- TRUE

dataset <- "mean"

my_metadata <- metadata[[dataset]]

target_data <- which(sapply(otu_table, nrow) <= 24)

# grab diversity data
temp <- do.call(
  "cbind.data.frame",  
  lapply(
    div_metrics, 
    function(x) {
      cbind.data.frame(x[target_data])
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

# set values for which we'd like to generate predictions
x_pred <- seq(from = 0, to = 5000, length = 101)

data_names <- colnames(my_data)[5:ncol(my_data)]

div_dist <- list()

for(i in 1:length(data_names)){

	div_dist[[data_names[i]]]$dat <- data.frame(
	  x = my_data$dist, 
	  logx = log(my_data$dist + 1), 
	  y = my_data[,data_names[i]]
	)

	div_dist[[data_names[i]]]$out <- lm(
	  y ~ x, 
      data = div_dist[[data_names[i]]]$dat
	)

	div_dist[[data_names[i]]]$conf <- predict(
	    object   = div_dist[[data_names[i]]]$out, 
	    newdata  = data.frame(x = x_pred), 
	    interval = "confidence", 
	    level    = 0.95
	)

	div_dist[[data_names[i]]]$outlog <- lm(
	  y ~ logx, 
      data = div_dist[[data_names[i]]]$dat
	)

	div_dist[[data_names[i]]]$conflog <- predict(
	    object   = div_dist[[data_names[i]]]$outlog, 
	    newdata  = data.frame(logx = x_pred), 
	    interval = "confidence", 
	    level    = 0.95
	)

}

# plot diversity as a function of distance from shore (full)

plot_name <- "diversity_distance_all"

if(!exists("legend_text")){ legend_text <- list()}
legend_text[plot_name] <- {
"Aggregate diversity metrics of each site plotted against distance from shore.
Both Simpson's Index (top) and richness (bottom) are shown for a variety of data subsets and transformations (left to right: mean, unfiltered mean, log(x + 1), transformed, scaled, spatially variable, and taxon clustered). 
Lines and bands illustrate the fit and 95% confidence interval of a linear model.
See methods text for detailed data descriptions."

}
if(EXPORT){
  pdf_file    <- file.path(fig_dir, paste(plot_name, ".pdf", sep = ""))
  legend_file <- file.path(fig_dir, paste(plot_name, "_legend.txt", sep = ""))
  writeLines(legend_text[[plot_name]], con = legend_file)
  pdf(file = pdf_file, width = 16, height = 8*(2/3))
}
par(mfrow = c(length(div_metrics), length(target_data)))
for(i in 1:length(div_dist)){
	par(mar = c(4,4,1,1))
	plot(x = div_dist[[i]]$dat$x, y = div_dist[[i]]$dat$y,
	  # log = "x", 
	  xlab = "Distance from shore (meters)", 
	  ylab = names(div_dist)[i]
	)
	plot_model(div_dist[[i]]$conf, x_pred, line_type = 2)
	# abline(div_dist[[i]]$out, col = "red")
}
if(EXPORT){
  dev.off()
}
dataset <- grep(".mean$", names(div_dist))
(summaries <- lapply(div_dist[dataset], summary))

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# plot diversity as a function of distance from shore, 
# for only a single dataset ("mean")

col_sub <- grep(".mean$", names(div_dist))

plot_name <- "diversity_distance"

if(!exists("legend_text")){ legend_text <- list()}
legend_text[plot_name] <- {
"Aggregate diversity metrics of each site plotted against distance from shore.
Both Simpson's Index (left) and richness (right) are shown, and have been computed from the mean abundance of unique DNA sequences found across 4 PCR replicates at each of 24 sites.
Lines and bands illustrate the fit and 95% confidence interval of a linear model.
"
}
if(EXPORT){
  pdf_file    <- file.path(fig_dir, paste(plot_name, ".pdf", sep = ""))
  legend_file <- file.path(fig_dir, paste(plot_name, "_legend.txt", sep = ""))
  writeLines(legend_text[[plot_name]], con = legend_file)
  pdf(file = pdf_file, width = 10, height = 5)
}

par(mfrow = c(1, length(col_sub)))

for(i in 1:length(col_sub)){

	par(mar = c(4,4,1,1))

	metric_name <- strsplit(names(div_dist)[col_sub[i]], "\\.")[[1]][1]
	
	# ylims <- range(c(range(div_dist[[metric_name]]$conf), range(temp_dat$y)))

	plot(x = div_dist[[col_sub[i]]]$dat$x, y = div_dist[[col_sub[i]]]$dat$y,
	  # log = "x", 
	  # axes = FALSE,
	  xlab = "Distance from shore (meters)", 
	  ylab = metric_name
	)

	# tick_x <- c(0,50,250,500,1000,2000,4000)
	# axis(1, at = log(tick_x + 1), labels = tick_x)
	
	# axis(2)
	
	# box()

	plot_model(div_dist[[col_sub[i]]]$conf, x_pred, line_type = 2)

}
if(EXPORT){
  dev.off()
}
