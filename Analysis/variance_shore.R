#!/usr/bin/env Rscript
#-------------------------------------------------------------------------------

EXPORT <- FALSE
# choose which dataset to use
# from names(otu_table), (mean, spvar, mean_unfilt, scale01)
dataset     <- "mean"
use_binary  <- FALSE
my_table    <- otu_table[[dataset]]
my_metadata <- metadata[[dataset]]


#-------------------------------------------------------------------------------
data_by_dist <- split(
	x = as.data.frame(my_table),
	f = my_metadata[,colname_ycoord]
)

distances <- rep(names(data_by_dist), times = sapply(data_by_dist, nrow))

dis_by_dist <- lapply(data_by_dist, function(x)
	1 - as.vector(vegdist(x, method = "bray", binary = use_binary))
)

my_data <- data.frame(x = as.numeric(distances), y = unlist(dis_by_dist))

# transform distance to shore

# run linear model
sim_from_shore <- list()

sim_from_shore[["linear"]]$out <- lm(y ~ log(x+1), data = my_data)

# no significant slope
summary(sim_from_shore[["linear"]]$out)


sim_from_shore[["MM_asy0"]] <- list()

(sim_from_shore[["MM_asy0"]]$out <- nls(
  formula   = models[["MM_asy0"]]$form, # from distance_decay.R
  data      = my_data,
  start     = models[["MM_asy0"]]$init, # from distance_decay.R
  lower     = params["min",],
  upper     = params["max",],
  algorithm = "port"
))
summary(sim_from_shore[["MM_asy0"]]$out)

x_pred <- seq(from = 0, to = 5000, length = 101)

temp <- predictNLS(
  model    = sim_from_shore[["MM_asy0"]]$out,
  newdata  = data.frame(x = x_pred),
  interval = "confidence",
  alpha    = 0.05
)$summary # the output of predictNLS is huge, just keep the summary
sim_from_shore[["MM_asy0"]]$conf <- data.frame(
  fit = temp[,"Prop.Mean.1"],
  lwr = temp[,5], # "Prop.2.5%", but varies depending on alpha level
  upr = temp[,6]  # "Prop.97.5%"
)
rm(temp)
sim_from_shore[["MM_asy0"]]$pred <- sim_from_shore[["MM_asy0"]]$conf[,"fit"]

#-------------------------------------------------------------------------------
# PLOTTING
# initiate a list to store legend text if it doesn't already exist
if(!exists("legend_text")){
  legend_text <- list()
}

plot_name   <- "sim_shore"

legend_text[plot_name] <- {
"Similarity of environmental DNA communities by distance from and along shore.
Each point represents the Bray-Curtis similarity among sites at the same distance along three parallel transects comprising a 3000 by 4000 meter grid, and red diamonds represent within-distance means.
A linear model of the log-transformed distance could not distinguish from a slope of 0."
}

if(EXPORT){
  pdf_file    <- file.path(fig_dir, paste(plot_name, ".pdf", sep = ""))
  legend_file <- file.path(fig_dir, paste(plot_name, "_legend.txt", sep = ""))
  writeLines(legend_text[[plot_name]], con = legend_file)
  pdf(file = pdf_file, width = 8, height = 4)
}
par(mfrow = c(1,2))
par(mar = c(4,4,1,1))
plot(
	x = my_data$x,
	xlab = "Distance from shore",
	y = my_data$y,
	ylab = "Community Similarity",
	ylim = c(0,1),
	pch = 21,
	col = hsv(1,1,0),
	bg  = hsv(1,1,0,0.2),
	log = "x",
	las = 1
)
plot_model(sim_from_shore[["MM_asy0"]]$conf, pred_vec = x_pred)
points(
	y = sapply(dis_by_dist, mean),
	x = as.numeric(names(dis_by_dist)),
	pch = 23,
	col = hsv(1,1,1),
	bg  = hsv(1,1,1,0.2)
)
# abline(lm_out, col = hsv(0.6, 1, 1), lwd = 2, lty = 2)

#-------------------------------------------------------------------------------
# ALONG SHORE
data_by_dist <- split(
	x = as.data.frame(my_table),
	f = my_metadata[,colname_xcoord]
)

dis_by_dist <- lapply(data_by_dist, function(x)
	1 - as.vector(vegdist(x, method = "bray", binary = use_binary))
)

distances <- rep(names(dis_by_dist), times = sapply(dis_by_dist, length))

my_data <- data.frame(
  x = as.numeric(distances), #, times = sapply(dis_by_dist, length))
  y = unlist(dis_by_dist)
)

# transform distance to shore

# run linear model
sim_along_shore <- list()

sim_along_shore[["linear"]]$out <- lm(y ~ log(x+1), data = my_data)

# no significant slope
summary(sim_along_shore[["linear"]]$out)


sim_along_shore[["MM_asy0"]] <- list()

(sim_along_shore[["MM_asy0"]]$out <- nls(
  formula   = models[["MM_asy0"]]$form, # from distance_decay.R
  data      = my_data,
  start     = models[["MM_asy0"]]$init, # from distance_decay.R
  lower     = params["min",],
  upper     = params["max",],
  algorithm = "port"
))
summary(sim_from_shore[["MM_asy0"]]$out)

x_pred <- seq(from = 0, to = 2000, length = 10)

temp <- predictNLS(
  model    = sim_along_shore[["MM_asy0"]]$out,
  newdata  = data.frame(x = x_pred),
  interval = "confidence",
  alpha    = 0.05
)$summary # the output of predictNLS is huge, just keep the summary
sim_along_shore[["MM_asy0"]]$conf <- data.frame(
  fit = temp[,"Prop.Mean.1"],
  lwr = temp[,5], # "Prop.2.5%", but varies depending on alpha level
  upr = temp[,6]  # "Prop.97.5%"
)
rm(temp)
sim_along_shore[["MM_asy0"]]$pred <- sim_along_shore[["MM_asy0"]]$conf[,"fit"]

par(mar = c(4,4,1,1))
plot(
	x = my_data$x,
	xlab = "Distance along shore",
	y = my_data$y,
	ylab = "Community Similarity",
	ylim = c(0,1),
	pch = 21,
	col = hsv(1,1,0),
	bg  = hsv(1,1,0,0.2),
	# log = "x",
	las = 1
)
plot_model(sim_along_shore[["MM_asy0"]]$conf, pred_vec = x_pred)
points(
	y = sapply(dis_by_dist, mean),
	x = as.numeric(names(dis_by_dist)),
	pch = 23,
	col = hsv(1,1,1),
	bg  = hsv(1,1,1,0.2)
)
if(EXPORT){
  dev.off()
}
# END ALONG SHORE
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# individual otu abundance variance
var_by_dist <- lapply(data_by_dist, function(x)
	apply(x, MARGIN = 2, FUN = function(x) var(log(x + 1)))
)

distances <- rep(names(var_by_dist), times = sapply(var_by_dist, length))

to_plot <- data.frame(dist = as.numeric(distances), var = unlist(var_by_dist))

#-------------------------------------------------------------------------------
# PLOTTING
# initiate a list to store legend text if it doesn't already exist
if(!exists("legend_text")){
  legend_text <- list()
}

plot_name   <- "variance_from_shore"

legend_text[plot_name] <- {
"Variance among transects with distance from shore for each OTU."
}

if(EXPORT){
  pdf_file    <- file.path(fig_dir, paste(plot_name, ".pdf", sep = ""))
  legend_file <- file.path(fig_dir, paste(plot_name, "_legend.txt", sep = ""))
  writeLines(legend_text[[plot_name]], con = legend_file)
  pdf(file = pdf_file, width = 4, height = 4)
}
par(mar = c(4,4,1,1))
plot(
	x = to_plot$dist,
	xlab = "Distance from shore",
	y = to_plot$var,
	ylab = "Within-OTU Variance of counts",
	# ylim = c(0,1),
	log = "x",
	pch = 21,
	col = hsv(1,1,0),
	bg  = hsv(1,1,0,0.2),
	las = 1
)
points(
	y = sapply(var_by_dist, mean),
	x = as.numeric(names(var_by_dist)),
	pch = 23,
	col = hsv(1,1,1),
	bg  = hsv(1,1,1,0.2)
)

# abline(lm_out, col = hsv(0.6, 1, 1), lwd = 2, lty = 2)
# pdf(file = file.path(fig_dir, "variance_from_shore.pdf"))
par(mar = c(4,4,1,1))
stripchart(
	x = var_by_dist,
	ylab = "Within-OTU Variance of counts",
	xlab = "Distance from shore",
	method = "jitter",
	pch = 21,
	col = hsv(1,1,0),
	bg  = hsv(1,1,0,0.2),
	las = 1,
	vertical = TRUE
)
points(
	x = sapply(var_by_dist, mean),
	pch = 23,
	col = hsv(1,1,1),
	bg  = hsv(1,1,1,0.2)
)
if(EXPORT){
  dev.off()
}
