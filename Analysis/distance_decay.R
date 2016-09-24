# FIRST RUN 'load_data.r'

# dissimilarity by distance for reduced data set (mean abundance within sample for each OTU)
#-------------------------------------------------------------------------------

# REQUIRES:
# 1. dataframe called "metadata" with unique sequenced samples given in column "sample_id"
my_metadata <- metadata[["mean"]] # metadata[!duplicated(metadata[,colname_env_sample]),]
# 2. OTU table with rownames that correspond to aforementioned column "sample_id"
my_table <- otu_table[["mean"]] # mean, mean_unfilt, spvar, otu_named, as.binary(otu_mean), log, filt
rownames(my_table) # should be e.g. PCT-C-0000 etc, aka "env_sample_name"

EXPORT <- FALSE # export plots/files?

library(geosphere) # distm()
library(vegan) # vegdist()
library(propagate) # predictNLS

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

# calculate replicate PCR similarity
otu_temp <- as.data.frame(prop(otu_table[["clean"]]))
PCR_similarities <- lapply(
  split(otu_temp, metadata[["clean"]][ , colname_env_sample]),
  function(x) {1 - vegdist(x, method = vegdist_method)}
)
rm(otu_temp)

#-------------------------------------------------------------------------------
# compare similarities of PCR replicates to environmental samples
PCR_mean <- round(mean(unlist(PCR_similarities)), digits = 3)
PCR_sd   <- round(  sd(unlist(PCR_similarities)), digits = 3)
env_mean <- round(mean(unlist(comm_dist["Bray_Curtis"])), digits = 3)
env_sd   <- round(  sd(unlist(comm_dist["Bray_Curtis"])), digits = 3)

text_similarity <- paste(
"PCR replicates within an environmental sample were extremely similar (",
PCR_mean, " plusminus ", PCR_sd,
") and far more similar than environmental samples (",
env_mean, " plusminus ", env_sd, ").",
sep = "")
par(mar = c(2,4,1,1))
boxplot(
list(PCR = unlist(PCR_similarities), environment = unlist(comm_dist["Bray_Curtis"])), 
las = 1, ylab = "Similarity"
)


#-------------------------------------------------------------------------------
# arrange the data for model fitting
model_data_full <- data.frame(
  dist2df(geo_dist),
  data.frame(lapply(comm_dist, as.vector))
)
# order the columns
model_data_full <- model_data_full[order(model_data_full[,"dist"]), ]

if(EXPORT){
# write the data
  write.csv(
    x = model_data_full,
    file = file.path(data_dir, "model_data_full.csv"),
    row.names = FALSE
  )
}

# subset for a given run of the models
model_data <- data.frame(
  x = model_data_full[,"dist"],
  y = model_data_full[,distance_name]
)


#-------------------------------------------------------------------------------
# Fit some models, estimate some parameters

# start a data frame to store parameter values,
# including initialization values (i.e. good guesses)
params <- data.frame(

# intercept
intercept = c(
  min  = 0,
  max  = 1,
  init = 0.5
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
  init = 0
),

# halflife
halflife = c(
  min  = 0,
  max  = Inf,
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

# This statement was under the linear model; can't remember where I got it
# "Then the regression coefficient is usually used in the literature as the descriptor of distance decay, or the distance at which 50% of the maximum similarity is observed."

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
    params["init",c("intercept", "asymptote", "halflife")]
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
    params["init",c("asymptote", "halflife")]
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
    params["init", c("intercept", "halflife")]
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
    c(halflife = params["init", "halflife"]) # required to get name
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
  init =
    params["init",c("intercept", "deltay", "halflife")]
)

# distinguish linear from nonlinear
which_linear    <- c("linear", "loglinear")
which_nonlinear <- names(models)[!names(models) %in% which_linear]

#-------------------------------------------------------------------------------
# run the linear models
for(i in which_linear){
  models[[i]]$out <- lm(formula = models[[i]]$form, data = model_data)
}
#-------------------------------------------------------------------------------
# run the nonlinear models
for(i in which_nonlinear){
  models[[i]]$out <- nls(
    formula      = models[[i]]$form,
    start        = models[[i]]$init,
    data         = model_data,
    lower        = params["min",],
    upper        = params["max",],
    algorithm    = "port"
  )
}

#-------------------------------------------------------------------------------
# summarize the models
# lapply(lapply(models, "[[", "out"), summary)
for(i in 1:length(models)){
  models[[i]]$summary <- summary(models[[i]]$out)
}

#-------------------------------------------------------------------------------
# run the predictions
x_pred <- seq(from = 0, to = 5000, length = 101)
for(i in 1:length(models)){
  models[[i]]$pred <- predict(models[[i]]$out,
    newdata = data.frame(x = x_pred))
}

# add confidence intervals to the predictions
# for the linear models
for(i in which_linear){
  models[[i]]$conf <- predict(
    object   = models[[i]]$out, 
    newdata  = data.frame(x = x_pred), 
    interval = "confidence", 
    level    = 0.95
  )
}
#-------------------------------------------------------------------------------
# and for the nonlinear models
# I can't get predictNLS to work for nonlinear model with single parameter.
which_pred <- which_nonlinear[which_nonlinear != "MM_int1asy0"]
for(i in which_pred){
  temp <- predictNLS(
    model    = models[[i]]$out,
    newdata  = data.frame(x = x_pred),
    interval = "confidence",
    alpha    = 0.05
  )$summary # the output of predictNLS is huge, just keep the summary
  models[[i]]$conf <- data.frame(
    fit = temp[,"Prop.Mean.1"], 
    lwr = temp[,5], # "Prop.2.5%", but varies depending on alpha level
    upr = temp[,6]  # "Prop.97.5%" 
  )
  rm(temp)
}



#-------------------------------------------------------------------------------
# Other analyses to consider:

# Mantel Test
# mantel(comm_dist, geo_dist, perm = 9999)

# Generalized Additive Model
# library(mgcv) #gam
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# SAVE MODEL OUTPUT
if(EXPORT){
writeLines(
  capture.output(lapply(models, "[", c("out", "summary"))),
  con =  "model_output.txt"
)
}
# ALT:
# sink("lm.txt")
# summary(lm(cars$speed ~ cars$dist))
# sink()
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# PLOTTING
geo_dist_scaled <- log(model_data$x + 100)
plot_x <- model_data$x # geo_dist_scaled

if(EXPORT){
  plot_base   <- "distance_decay"
  pdf_file    <- file.path(fig_dir, paste(plot_base, ".pdf", sep = ""))
  legend_file <- file.path(fig_dir, paste(plot_base, "_legend.txt", sep = ""))
  writeLines(
"Distance decay relationship of environmental DNA communities. Each point represents the Bray-Curtis similarity of a site sampled along three parallel transects comprising a 3000 by 4000 meter grid. Blue dashed line represents fit of a nonlinear least squares regression (see Methods), and shading denotes the 95% confidence interval. Boxplot is comparisons within-sample across PCR replicates, separated by a vertical line at zero, where the central line is the median, the box encompasses the interquartile range, and the lines extend to 1.5 times the interquartile range. Boxplot outliers are omitted for clarity.",
  con = legend_file)
  pdf(file = pdf_file, width = 8, height = 4) #, width = 8, height = 3

}
par(mar = c(4,4,1,1))
plot(
	x = model_data$x,
	y = model_data$y,
	ylim = c(0,1),
	xlim = c(-500, max(model_data$x)),
	xaxt = "n",
	pch = 21,
	cex = 1,
	col = hsv(h = 0, s = 1, v = 0, alpha = 0.5),
	bg = rgb(0,0,0,alpha = 0.1 ), #,alpha = 0.1
	xlab = "Distance between samples (meters)",
	ylab = "Community similarity", #paste(,"(", distance_name, ")", sep = ""),
	# log = "x",
	axes = FALSE,
	las = 1
)
axis(side = 1, 
  at      = c(-350, 50, 1000, 2000, 3000, 4000), 
  labels = c("PCR", 50, 1000, 2000, 3000, 4000), lwd = 0, lwd.ticks = 1)
#, at = unique(log(metadata$dist_from_shore + 100)), labels = unique(metadata$dist_from_shore)
# abline(v = unique(log(metadata$dist_from_shore + 100)))
axis(side = 2, lwd = 0, lwd.ticks = 1, las = 1)
box()

# add PCR similarity
abline(v = 0)
boxplot(unlist(PCR_similarities), add = TRUE,
  at = -350,
  boxwex = 800,
  axes = FALSE,
  outpch = NA, # suppress plotting of outliers
  boxlwd = 0.5,
  medlwd = 0.5,
  show.names = FALSE
)

main_models <- c("MM_asy0") 

# add confidence band
x_bounds <- c(x_pred, rev(x_pred))
conf     <- models[[main_models]]$conf
y_bounds <- c(conf[,"lwr"], rev(conf[,"upr"]))
polygon(
  x = x_bounds, 
  y = y_bounds,
  col = hsv(h = 1, s = 1, v = 0.1, alpha = 0.2), 
  border = NA
)

# "linear","loglinear","MM_full","MM_int1","MM_asy0","MM_int1asy0","Harold" 
line_colors <- c("#6495ED") #6495ED, #0084d1 , "purple", "red"
line_types <- c(2)
for(i in 1:length(main_models)) {
	lines(x = x_pred, y = models[[main_models[[i]]]]$conf[,"fit"], 
      col = line_colors[i], lwd = 3, lty = line_types[i])
}

# legend(
  # "topright",
  # legend = names(model_pred)[which_models],
  # bty = "n", lty = line_types, col = line_colors, lwd = 3)

# title(main = "abundance", line = 0.1, adj = 0)

# Add LOESS line
# smoothed <- loess.smooth(
				# x = model_data$x,
				# y = model_data$y,
				# span = 2/3,
				# degree = 1,
				# family = "gaussian"
				# )
# lines(smoothed, col="red", lwd=2)
if(EXPORT){
  dev.off()
}

#===============================================================================
# plot all models
if(EXPORT){
  plot_pdf    <- paste(plot_base, "_all", ".pdf", sep = "")
  plot_pdf   <- file.path(fig_dir, plot_pdf)
  pdf(file = plot_pdf, width = 5, height = 5) #, width = 8, height = 3
}
for(i in 1:length(models)){
  par(mar = c(4,5,1,1))
  plot(
    x = model_data$x,
    y = model_data$y,
    ylim = c(0,1),
    xlim = c(0, max(model_data$x)),
    xaxt = "n",
    pch = 21,
    cex = 1,
    col = hsv(h = 0, s = 1, v = 0, alpha = 0.5),
    bg = rgb(0,0,0,alpha = 0.1 ), #,alpha = 0.1
    xlab = "Distance between samples (meters)",
	ylab = "Community similarity", #paste(,"(", distance_name, ")", sep = ""),
    # log = "x",
    axes = FALSE,
    las = 1
  )
  axis(side = 1, lwd = 0, lwd.ticks = 1)
  axis(side = 2, lwd = 0, lwd.ticks = 1, las = 1)
  box()

# add confidence band
  if(length(models[[i]]$conf) > 0){
    x_bounds <- c(x_pred, rev(x_pred))
    conf     <- models[[i]]$conf
    y_bounds <- c(conf[,"lwr"], rev(conf[,"upr"]))
    polygon(
      x = x_bounds, 
      y = y_bounds,
      col = hsv(h = 1, s = 1, v = 0.1, alpha = 0.2), 
      border = NA
    )
  }

  # add fit line
  lines(x = x_pred, y = models[[i]]$pred, col = "indianred", lwd = 3, lty = 3)
  legend("topright", legend = names(models)[i],
    bty = "n", lty = 3, col = "indianred", lwd = 3)
}
if(EXPORT){
  dev.off()
}
