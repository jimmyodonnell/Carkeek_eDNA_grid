#===============================================================================
## run distance decay analysis using a bunch of different datasets and metrics
## and plot multiple panels on a single file

EXPORT <- FALSE

data_set <- c("mean", "mean_unfilt", "log", "scale01", "spvar")
many_models <- list()


#-------------------------------------------------------------------------------
## Begin data subset loop
for(dat in 1:length(data_set)){

many_models[[data_set[dat]]] <- list()

# from names(otu_table)
# REQUIRES:
# 1. dataframe called "metadata" with unique sequenced samples given in column "sample_id"
my_metadata <- metadata[[data_set[dat]]]
# 2. OTU table with rownames that correspond to aforementioned column "sample_id"
my_table <- otu_table[[data_set[dat]]]

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
# calculate pairwise similarity of ecological communities

vegdist_methods <- c(
  "Bray_Curtis"   = "bray",
  "Morisita"      = "morisita",
  "Morisita_Horn" = "horn",
  "Jaccard"       = "jaccard",
  "Gower"         = "gower")

comm_sim <- list()
for(i in 1:length(vegdist_methods)){
  comm_sim[[names(vegdist_methods[i])]] <- vegdist(
    x = my_table, method = vegdist_methods[i], binary = FALSE)
}
for(i in 1:length(vegdist_methods)){
  comm_sim[[paste(names(vegdist_methods[i]), "bin", sep = "_")]] <- vegdist(
    x = my_table, method = vegdist_methods[i], binary = TRUE)
}

comm_sim <- lapply(comm_sim, function(x) 1 - x)

if(!(identical(attr(comm_sim[[1]], "Labels"), attr(geo_dist, "Labels")))){
	warning("Whoa there! the row/column names of the two distance matrices do not seem to add up. This is bad.")
}

#-------------------------------------------------------------------------------
# arrange the data for model fitting
model_data_full <- data.frame(
  dist2df(geo_dist),
  data.frame(lapply(comm_sim, as.vector))
)
# order the columns
model_data_full <- model_data_full[order(model_data_full[,"dist"]), ]

# if(EXPORT){
# # write the data
#   write.csv(
#     x = model_data_full,
#     file = file.path(data_dir, "model_data_full.csv"),
#     row.names = FALSE
#   )
# }

dist_met <- c("Bray_Curtis", "Bray_Curtis_bin", "Morisita_Horn", "Morisita_Horn_bin") # names(model_data_full)

#-------------------------------------------------------------------------------
## Begin similarity metric loop
for(metric in 1:length(dist_met)){

  # subset for a given run of the models
  model_data <- data.frame(
    x = model_data_full[,"dist"],
    y = model_data_full[,dist_met[metric]]
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
## bind to big list
many_models[[data_set[dat]]][[dist_met[metric]]] <- models

  } # for(metric in 1:length(dist_met))
  
  many_models[[data_set[dat]]]$model_data <- model_data_full

} # for(dat in 1:length(data_set)){


#-------------------------------------------------------------------------------
# Other analyses to consider:

# Mantel Test
# mantel(comm_sim, geo_dist, perm = 9999)

# Generalized Additive Model
# library(mgcv) #gam
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# SAVE MODEL OUTPUT
# if(EXPORT){
# writeLines(
#   capture.output(lapply(models, "[", c("out", "summary"))),
#   con =  "model_output.txt"
# )
# }
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# PLOTTING

plot_name   <- "distance_decay_multi"

if(!exists("legend_text")){
  legend_text <- list()
}


legend_text[plot_name] <- "
Distance decay relationship of environmental DNA communities using a variety of models, metrics, and data subsets.
Each point represents the similarity of a site sampled along three parallel transects comprising a 3000 by 4000 meter grid.
Rows of plots represent different data subsets, while columns represent different similarity metrics.
Results using the Jaccard distance are omitted because of its similarity to Bray-Curtis."

if(EXPORT){
  pdf_file    <- file.path(fig_dir, paste(plot_name, ".pdf", sep = ""))
  legend_file <- file.path(fig_dir, paste(plot_name, "_legend.txt", sep = ""))
  writeLines(legend_text[[plot_name]], con = legend_file)
  pdf(file = pdf_file, width = 8, height = 10)
}

subset_names <- c("Filtered", "Unfiltered", "Log", "Scaled {0,1}", 
  "Spatially Variable")

metric_names <- c("Bray-Curtis", "Bray-Curtis (binary)", 
  "Morisita-Horn", "Morisita-Horn (binary)")

par(
  oma = c(2,2,2,2), # set outer margins
  mfrow = c(length(many_models), length(dist_met)) # set dimensions
)
for(subset in 1:length(many_models)){
  for(metric in 1:length(dist_met)){
    plot_points_nolab(
      x = many_models[[subset]]$model_data[,"dist"], 
      y = many_models[[subset]]$model_data[,dist_met[metric]]
    )
    plot_model(
      model_list_item = many_models[[subset]][[dist_met[metric]]][["MM_asy0"]], 
      pred_vec = x_pred, 
      line_color = hsv(1,1,1), 
      band_color = hsv(1,1,1, alpha = 0.3) 
    )
    if(subset == 1){
      mtext(metric_names[metric], side = 3, line = 1)
    }
    if(metric == length(dist_met)){
      mtext(subset_names[subset], side = 4, line = 1)
    }
  }
}

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(4, 4, 1, 4), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", ann = FALSE)
mtext("Community Similarity", side = 2, line = 2, cex = 1.5)
mtext("Distance between samples (meters) ", side = 1, line = 3, cex = 1.5)


# legend(
  # "topright",
  # legend = names(model_pred)[which_models],
  # bty = "n", lty = line_types, col = line_colors, lwd = 3)

if(EXPORT){
  dev.off()
}
