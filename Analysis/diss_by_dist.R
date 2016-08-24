# FIRST RUN 'load_data.r'

# dissimilarity by distance for reduced data set (mean abundance within sample for each OTU)
#-------------------------------------------------------------------------------

# REQUIRES:
# 1. dataframe called "metadata" with unique sequenced samples given in column "sample_id"
my_metadata <- metadata_mean # metadata[!duplicated(metadata[,colname_env_sample]),]
# 2. OTU table with rownames that correspond to aforementioned column "sample_id"
my_table <- otu_filt # otu_mean, otu_spvar, otu_named, as.binary(otu_mean), otu_log, otu_filt
rownames(my_table) # should be e.g. PCT-C-0000 etc, aka "env_sample_name"


library(geosphere) # distm()

export_plots <- FALSE

vegdist_method <- "bray"

distance_name <- switch(vegdist_method,
       bray     = "Bray-Curtis", 
       morisita = "Morisita", 
       jaccard  = "Jaccard", 
       gower    = "Gower")

# calculate pairwise great circle distance between sampling locations using Haversine method
geo_dist <- as.dist(distm(x = my_metadata[,c(colname_lon, colname_lat)], fun = distHaversine))
geo_dist_v <- as.vector(geo_dist)
# dimnames(geo_dist) <- list(my_metadata$env_sample_name, my_metadata$env_sample_name)

# vegdist_bin <- c(TRUE, FALSE)
# for(i in 1:length(vegdist_bin)){
comm_dist <- vegdist(my_table, method = vegdist_method, binary = FALSE) #, diag = TRUE, upper = TRUE
comm_dist_v <- as.vector(comm_dist)

if(!(identical(dimnames(comm_dist), dimnames(geo_dist)))){
	warning("Whoa there! the row/column names of the two distance matrices do not seem to add up. This is bad.")
}

# arrange the data for model fitting
the_data <- data.frame(comm = as.vector(comm_dist), space = as.vector(geo_dist))

#-------------------------------------------------------------------------------
# Fit some models, estimate some parameters
model_out <- list()
model_pred <- list()
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Michaelis-Menten
fix_asymptote_1 <- FALSE
if(fix_asymptote_1){
	Vm <- 1
	start_list <- list(Km = max(the_data$comm)/2)
} else {
	start_list <- list(Vm = max(the_data$com), Km = max(the_data$comm)/2)
}
mm_fit <- nls(
  formula = comm ~ Vm * space/(Km + space),
  data = the_data,
  start = start_list)
pred_mm <- predict(mm_fit)
model_out[["Michaelis-Menten"]] <- mm_fit
model_pred[["Michaelis-Menten"]] <- data.frame(x = sort(geo_dist), y = sort(pred_mm))

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Nonlinear Least Squares Regression (2 parameters)
nls_p2 <- nls(
	formula = comm_dist_v ~ 1 - ((1-INT) * exp( -RATE * geo_dist_v )),
	start = c(RATE = 0.02, INT = 0)
)
summary(nls_p2)
pred_nls_p2 <- predict(nls_p2)
model_out[["NLS-2p"]] <- nls_p2
model_pred[["NLS-2p"]] <- data.frame(x = sort(geo_dist), y = sort(pred_nls_p2))
# lines(sort(geo_dist), sort(pred_nls_p2), col = "purple", lty = 3, lwd = 2)
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Nonlinear Least Squares Regression (3 parameters)
nls_p3 <- nls(
	formula = comm_dist_v ~ INT + ASY_DIFF * (1 - exp( -RATE * geo_dist_v )),
	start = c(ASY_DIFF = 1, RATE = 0.02, INT = 0)
)
summary(nls_p3)
pred_nls_p3 <- predict(nls_p3)
model_out[["NLS-3p"]] <- nls_p3
model_pred[["NLS-3p"]] <- data.frame(x = sort(geo_dist), y = sort(pred_nls_p3))
# lines(sort(geo_dist), sort(pred_nls_p3), col = "indianred")
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Linear Model
lm_out <- lm(comm_dist ~ geo_dist)
lm_out <- lm(log(comm_dist)~ geo_dist)
# Then the regression coefficient is usually used in the literature as the descriptor of distance decay, or the distance at which 50% of the maximum similarity is observed.
summary(lm_out)
pred_lm <- predict(lm_out)
model_out[["Linear"]] <- lm_out
model_pred[["Linear"]] <- data.frame(x = sort(geo_dist), y = sort(pred_lm))
# lines(sort(geo_dist), sort(pred_lm), col = "blue")
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Other analyses to consider:

# Mantel Test
# mantel(comm_dist, geo_dist, perm = 9999)

# Generalized Additive Model
# library(mgcv) #gam 
# ?gam
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# SAVE MODEL OUTPUT
for(i in 1:length(model_out)){
  writeLines(
    capture.output(
      model_out[[i]], 
      summary(model_out[[i]])
      ), 
  con =  paste("model_output_", names(model_out[i]), ".txt", sep = ""))
}

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# PLOTTING
geo_dist_scaled <- log(geo_dist + 100)
plot_x <- geo_dist # geo_dist_scaled

if(export_plots){
  plot_base   <- "diss_by_dist"
  plot_pdf    <- paste(plot_base, ".pdf", sep = "")
  legend_file <- paste(plot_base, "_legend.txt", sep = "")
  writeLines(
"Distance decay relationship of environmental DNA communities. Each point represents a site sampled along three parallel transects comprising a 3000 by 4000 meter grid. Blue dashed line represents fit to a Michaelis-Menten function. Purple dotted line represents nonlinear least squares regression (see Methods).", 
             con = file.path(fig_dir, legend_file))
  pdf(file = file.path(fig_dir, plot_pdf)) #, width = 8, height = 3
  
}
par(mar = c(4,5,1,1))
plot(
	x = plot_x,
	y = comm_dist,
	ylim = c(0,1),
	xaxt = "n",
	pch = 21,
	cex = 1,
	col = hsv(h = 0, s = 1, v = 0, alpha = 0.5),
	bg = rgb(0,0,0,alpha = 0.1 ), #,alpha = 0.1
	xlab = "Distance between samples (meters)",
	ylab = paste("Community dissimilarity (", distance_name, ")", sep = ""), 
	# log = "x", 
	axes = FALSE,
	las = 1
)
axis(side = 1, lwd = 0, lwd.ticks = 1)
#, at = unique(log(metadata$dist_from_shore + 100)), labels = unique(metadata$dist_from_shore)
# abline(v = unique(log(metadata$dist_from_shore + 100)))
axis(side = 2, lwd = 0, lwd.ticks = 1, las = 1)
	
line_colors <- c("#6495ED", "purple", "red") #6495ED, #0084d1
line_types <- c(2,3)
for(model in c(1, 2)) {
	lines(model_pred[[model]], col = line_colors[model], lwd = 2, lty = line_types[model])
}

legend("bottomright", legend = names(model_pred)[c(1,2,4)], bty = "n", lty = line_types, col = line_colors, lwd = 2)

# title(main = "abundance", line = 0.1, adj = 0)

# Add LOESS line
# smoothed <- loess.smooth(
				# x = plot_x,
				# y = comm_dist,
				# span = 2/3,
				# degree = 1,
				# family = "gaussian"
				# )
# lines(smoothed, col="red", lwd=2)
if(export_plots){
  dev.off()
}
# }
