# FIRST RUN 'load_data.r'

# dissimilarity by distance for reduced data set (mean abundance within sample for each OTU)
#-------------------------------------------------------------------------------

# REQUIRES:
# 1. dataframe called "metadata" with unique sequenced samples given in column "sample_id"
my_metadata <- metadata[["mean"]] # metadata[!duplicated(metadata[,colname_env_sample]),]
# 2. OTU table with rownames that correspond to aforementioned column "sample_id"
my_table <- otu_table[["mean"]] # mean, mean_unfilt, spvar, otu_named, as.binary(otu_mean), log, filt
rownames(my_table) # should be e.g. PCT-C-0000 etc, aka "env_sample_name"

export_plots <- TRUE

library(geosphere) # distm()

#-------------------------------------------------------------------------------
# calculate pairwise great circle distance between sampling locations using Haversine method
geo_dist <- as.dist(
  distm(
    x = my_metadata[,c(colname_lon, colname_lat)], 
    fun = distHaversine
  )
)
attr(geo_dist, "Labels") <- my_metadata[, colname_env_sample]
# dimnames(geo_dist) <- list(my_metadata$env_sample_name, my_metadata$env_sample_name)

#-------------------------------------------------------------------------------
# calculate pairwise dis/similarity of ecological communities
vegdist_method <- "bray"

distance_name <- switch(vegdist_method,
       bray     = "Bray_Curtis", 
       morisita = "Morisita", 
       horn     = "Morisita-Horn", 
       jaccard  = "Jaccard", 
       gower    = "Gower")

vegdist_methods <- c(
  "Bray_Curtis"   = "bray", 
  "Morisita"      = "morisita", 
  "Morisita_Horn" = "horn", 
  "Jaccard"       = "jaccard", 
  "Gower"         = "gower")

USE_SIMILARITY <- TRUE # use similarity instead of dissimilarity

# vegdist_bin <- c(TRUE, FALSE)
# for(i in 1:length(vegdist_bin)){

comm_dist <- list()
for(i in 1:length(vegdist_methods)){
  comm_dist[[names(vegdist_methods[i])]] <- vegdist(
    x = my_table, method = vegdist_methods[i], binary = FALSE)
}
for(i in 1:length(vegdist_methods)){
  comm_dist[[paste(names(vegdist_methods[i]), "bin", sep = "_")]] <- vegdist(
    x = my_table, method = vegdist_methods[i], binary = TRUE)
}

if(USE_SIMILARITY){
  comm_dist <- lapply(comm_dist, function(x) 1 - x)
}

if(!(identical(attr(comm_dist[[1]], "Labels"), attr(geo_dist, "Labels")))){
	warning("Whoa there! the row/column names of the two distance matrices do not seem to add up. This is bad.")
}

# arrange the data for model fitting
model_data_full <- data.frame(
  dist2df(geo_dist), 
  data.frame(lapply(comm_dist, as.vector))
)
# order the columns
model_data_full <- model_data_full[order(model_data_full[,"dist"]), ]
# write the data
write.csv(
  x = model_data_full, 
  file = file.path(data_dir, "model_data_full.csv"), 
  row.names = FALSE
)

# subset for a given run of the models
model_data <- data.frame(
  x = model_data_full[,"dist"], 
  y = model_data_full[,distance_name]
)


#-------------------------------------------------------------------------------
# Fit some models, estimate some parameters

# start a data frame to store parameter initialization values (i.e. good guesses)
params <- data.frame(

# intercept
intercept = c(
  min  = 0, 
  max  = 1, 
  init =   0.5
), 

# intercept_mod
int_mod = c( # I don't think I should use this parameter/formulation
  min  = -Inf, 
  max  = Inf, 
  init = NA
), 

# asymptote
asymptote = c( 
  min  = 0, 
  max  = 1, 
  init =   0
),

# halflife
halflife = c(
  min  = 0, 
  max  = max(model_data$y), 
  init = max(model_data$y)/2
),

# delta y
deltay = c(
  min  = 0, 
  max  = 1, 
  init = 1
),

# slope
slope = c(
  min  = -Inf, 
  max  = Inf, 
  init = -0.001
)

)

# start a vector to store parameter minimum values

#-------------------------------------------------------------------------------
# set up formula and equations for each model
models   <- list()

#===============================================================================
# Linear Model
models[["linear"]] <- list(
  func =
  function(x, intercept = 0, slope = 1){
    y <- intercept + slope * x
    return(y)
  },
  form =
    y  ~ x,
  init = NA
)

#===============================================================================
# Log-Linear Model
models[["loglinear"]] <- list(
  func =
  function(x, intercept = 0, slope = 1){
    y <- intercept + slope * log(x)
    return(y)
  },
  form =
    y  ~ log(x),
  init = NA
)

#===============================================================================
# Michaelis-Menten, full
models[["MM_full"]] <- list(
  func =
  function(x, intercept = 1, asymptote = 0, halflife = max(x)/2){
    int_mod <- intercept * halflife
    y <- ((intercept * halflife) + asymptote * x)/( halflife + x)
    return(y)
    },
  form =
    y  ~ ((intercept * halflife) + asymptote * x)/( halflife + x),
  init =
    init[c("intercept", "asymptote", "halflife")]
)

#===============================================================================
# Michaelis-Menten, intercept = 1
models[["MM_int1"]] <- list(
  func =
  function(x, asymptote = 0, halflife = max(x)/2){
    y <- (halflife + asymptote * x)/( halflife + x)
    return(y)
  },
  form =
    y  ~ (halflife + asymptote * x)/( halflife + x),
  init =
    init[c("asymptote", "halflife")]
)

#===============================================================================
# Michaelis-Menten, asymptote = 0
models[["MM_asy0"]] <- list(
  func =
  function(x, intercept = 1, halflife = max(x)/2){
    int_mod <- intercept * halflife
    y <- (intercept * halflife)/(halflife + x)
    return(y)
  },
  form =
    y  ~ (intercept * halflife)/(halflife + x),
  init =
    init[c("intercept", "halflife")]
)

#===============================================================================
# Michaelis-Menten, intercept = 1; asymptote = 0
models[["MM_int1asy0"]] <- list(
  func =
  function(x, halflife = max(x)/2){
    y <- halflife/(halflife + x)
    return(y)
  },
  form =
    y  ~ halflife/(halflife + x),
  init =
    init["halflife"]
)

#===============================================================================
# I have no idea what to call this weird looking model
models[["Harold"]] <- list(
  func =
    function(x, intercept = 1, deltay = -1, halflife = max(x)/2){
    y <- intercept + (deltay*x)/(halflife + x)
    return(y)
  },
  form =
    y  ~ intercept + (deltay*x)/(halflife + x),
  init = init[c("intercept", "deltay", "halflife")]
)

# distinguish linear from nonlinear
which_linear    <- c("linear", "loglinear")
which_nonlinear <- names(models)[!names(models) %in% which_linear]

#-------------------------------------------------------------------------------
# run the linear models
for(i in c("linear", "loglinear")){
  models[[i]]$out <- lm(formula = models[[i]]$form, data = model_data)
}
#-------------------------------------------------------------------------------
# run the nonlinear models
for(i in which_nonlinear){
  models[[i]]$out <- nls(
    formula      = models[[i]]$form, 
    start        = models[[i]]$init, 
    data         = model_data, 
    lower        = -Inf,
    upper        = Inf,
    algorithm    = "port"
  )
}

#-------------------------------------------------------------------------------
# summarize the models

lapply(lapply(models, "[[", "out"), summary)

#-------------------------------------------------------------------------------
# run the predictions
for(i in 1:length(models)){
  models[[i]]$pred <- predict(models[[i]]$out)
}

# add confidence intervals to the predictions

#-------------------------------------------------------------------------------
# Then the regression coefficient is usually used in the literature as the descriptor of distance decay, or the distance at which 50% of the maximum similarity is observed.
summary(lm_out)
model_out[["Linear"]] <- lm_out
#-------------------------------------------------------------------------------



model_out <- list()
model_pred <- list()
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Michaelis-Menten
if(!USE_SIMILARITY){
fix_asymptote <- TRUE
if(fix_asymptote){
	Vm <- 1
	start_list <- list(Km = max(comm_dist_v)/2)
} else {
	start_list <- list(Vm = max(comm_dist_v), Km = max(comm_dist_v)/2)
}
mm_fit <- nls(
  formula = comm_dist_v ~ Vm * geo_dist_v/(Km + geo_dist_v),
  start = start_list)
pred_mm <- predict(mm_fit)
model_out[["Michaelis-Menten"]] <- mm_fit
model_pred[["Michaelis-Menten"]] <- data.frame(
  x = c(0, sort(geo_dist_v)), 
  y = c(0, sort(pred_mm))
)
}
if(USE_SIMILARITY){
fix_asymptote <- FALSE
if(fix_asymptote){
	Vm <- 1
	start_list <- list(Km = max(comm_dist_v)/2)
} else {
	start_list <- list(Vm = max(comm_dist_v), Km = max(comm_dist_v)/2)
}
mm_fit <- nls(
  formula = comm_dist_v ~ 1 - (Vm * geo_dist_v)/(Km + geo_dist_v),
  start = start_list)
pred_mm <- predict(mm_fit)
model_out[["Michaelis-Menten"]] <- mm_fit
model_pred[["Michaelis-Menten"]] <- data.frame(
  x = c(0, sort(geo_dist)), 
  y = c(1, sort(pred_mm, decreasing = TRUE))
)
}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Nonlinear Least Squares Regression (2 parameters)
if(!USE_SIMILARITY){
nls_p2 <- nls(
	formula = comm_dist_v ~ 1 - ((1-INT) * exp( -RATE * geo_dist_v )),
	start = c(RATE = 0.02, INT = 0)
)
summary(nls_p2)
pred_nls_p2 <- predict(nls_p2)
model_out[["NLS-2p"]] <- nls_p2
model_pred[["NLS-2p"]] <- data.frame(x = sort(geo_dist), y = sort(pred_nls_p2))
}
# lines(sort(geo_dist), sort(pred_nls_p2), col = "purple", lty = 3, lwd = 2)
if(USE_SIMILARITY){
nls_p2 <- nls(
    formula = comm_dist_v ~ INT/(RATE + geo_dist_v),
    start   = c(INT = 1000, RATE = 1000)
	# formula = comm_dist_v ~ (1-INT) * exp( -RATE * geo_dist_v ),
	# start = c(RATE = 0.02, INT = 0)
)
summary(nls_p2)
pred_nls_p2 <- predict(nls_p2)
model_out[["NLS-2p"]] <- nls_p2
model_pred[["NLS-2p"]] <- data.frame(
  x = sort(geo_dist), 
  y = sort(pred_nls_p2, decreasing = TRUE)
)
}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Nonlinear Least Squares Regression (3 parameters)
if(!USE_SIMILARITY){
nls_p3 <- nls(
	formula = y ~ INT + ASY_DIFF * (1 - exp( -RATE * x )),
	data = model_data, 
	start = c(ASY_DIFF = 1, RATE = 0.02, INT = 0)
)
summary(nls_p3)
pred_nls_p3 <- predict(nls_p3)
model_out[["NLS-3p"]] <- nls_p3
model_pred[["NLS-3p"]] <- data.frame(x = sort(geo_dist), y = sort(pred_nls_p3))
# lines(sort(geo_dist), sort(pred_nls_p3), col = "indianred")
}
if(USE_SIMILARITY){
nls_p3 <- 
nls(
  formula = y ~ intercept / (halflife + x), 
  data = model_data, 
  # lower = c(INT = -Inf,   ASY_DIFF = 0, RATE = -Inf),
  # start   = c(intercept = 500, halflife = 1000), 
  start   = params["init",c("intercept", "halflife")], 
  # formula = comm_dist_v ~ INT - ASY_DIFF * (1 - exp(-RATE * geo_dist_v)), 
  # start = list(INT = 0.5, ASY_DIFF = 0.5, RATE = 0.02), 
  # lower = list(INT = 0, ASY_DIFF = 0.1,   RATE = -Inf), 
  # upper = list(INT = 1, ASY_DIFF = 0.5,   RATE = Inf), 
  algorithm = "port"
)
summary(nls_p3)
pred_nls_p3 <- predict(nls_p3)
model_out[["NLS-3p"]] <- nls_p3
model_pred[["NLS-3p"]] <- data.frame(x = sort(geo_dist), y = sort(pred_nls_p3, decreasing = TRUE))
}
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
# ALT:
# sink("lm.txt")
# summary(lm(cars$speed ~ cars$dist))
# sink()
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
"Distance decay relationship of environmental DNA communities. Each point represents a site sampled along three parallel transects comprising a 3000 by 4000 meter grid. Blue dashed line represents nonlinear least squares regression (see Methods). Boxplot is comparisons within-sample across PCR replicates, separated by a vertical line at zero, where the central line is the median, the box encompasses the interquartile range, and the lines extend to 1.5 times the interquartile range. Boxplot outliers are omitted for clarity.", 
             con = file.path(fig_dir, legend_file))
  pdf(file = file.path(fig_dir, plot_pdf), width = 8, height = 4) #, width = 8, height = 3
  
}
par(mar = c(4,5,1,1))
plot(
	x = plot_x,
	y = comm_dist,
	ylim = c(0,1),
	xlim = c(-500, max(plot_x)), 
	xaxt = "n",
	pch = 21,
	cex = 1,
	col = hsv(h = 0, s = 1, v = 0, alpha = 0.5),
	bg = rgb(0,0,0,alpha = 0.1 ), #,alpha = 0.1
	xlab = "Distance between samples (meters)",
	ylab = paste("Community similarity (", distance_name, ")", sep = ""), 
	# log = "x", 
	axes = FALSE,
	las = 1
)
axis(side = 1, at = c(-350, 50, 1000, 2000, 3000, 4000), labels = c("PCR", 50, 1000, 2000, 3000, 4000), lwd = 0, lwd.ticks = 1)
#, at = unique(log(metadata$dist_from_shore + 100)), labels = unique(metadata$dist_from_shore)
# abline(v = unique(log(metadata$dist_from_shore + 100)))
axis(side = 2, lwd = 0, lwd.ticks = 1, las = 1)
box()

# add PCR similarity
otu_temp <- as.data.frame(otu_clean/rowSums(otu_clean))
PCR_similarities <- lapply(
  split(otu_temp, metadata_clean[ , colname_env_sample]), 
  vegdist, method = vegdist_method
)
abline(v = 0)
boxplot(1 - unlist(PCR_similarities), add = TRUE, 
  at = -350, 
  boxwex = 800, 
  axes = FALSE,
  outpch = NA, # suppress plotting of outliers
  boxlwd = 0.5, 
  medlwd = 0.5, 
  show.names = FALSE
)

which_models <- c(1)
line_colors <- c("#6495ED", "purple", "red") #6495ED, #0084d1
line_types <- c(2)
for(model in which_models) {
	lines(model_pred[[model]], col = line_colors[model], lwd = 3, lty = line_types[model])
}

# legend(
  # "bottomright", 
  # legend = names(model_pred)[which_models], 
  # bty = "n", lty = line_types, col = line_colors, lwd = 3)

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
