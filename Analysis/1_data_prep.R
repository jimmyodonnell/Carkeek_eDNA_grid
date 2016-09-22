
EXPORT <- FALSE # export figures to files?

R_files <- list.files(path = "functions", pattern = "\\.R$", full.names = TRUE)
sapply(R_files, source)

# calculate proportional abundance of OTUs in each sample
# otu_table[["prop"]] <- otu_table[["raw"]]/rowSums(otu_table[["raw"]])
# metadata[["prop"]]  <- metadata[["raw"]]

#-------------------------------------------------------------------------------
CHECK_FOR_OUTLIERS <- TRUE
# If you'd like to check for and remove inconsistent PCR replicates, go to:
if(EXPORT){
  pdf(file = file.path(fig_dir, "PCR_consistency.pdf"))
}
if(CHECK_FOR_OUTLIERS) {
  cleaned <- find_bad_replicate(
  	my_table         = otu_table[["raw"]],
  	my_metadata      = metadata[["raw"]],
  	sample_id_column = colname_sampleid,
  	grouping_column  = colname_env_sample, 
  	threshold_sd     = 1, 
    plot_results     = TRUE
  )
  otu_table[["clean"]] <- strip_absent(cleaned[[1]])
  metadata[["clean"]]  <- cleaned[[2]]
  rm(cleaned)
}
if(EXPORT){
  dev.off()
}
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# RESCALE TO EQUAL SEQUENCING DEPTHS PER SAMPLE
RESCALE_SEQUENCING_DEPTH <- TRUE

if(RESCALE_SEQUENCING_DEPTH) {
	otu_table[["scaled"]] <- rescale_rowsums(otu_table[["clean"]])
    metadata[["scaled"]]  <- metadata[["clean"]]
}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# CALCULATE THRESHOLD AT WHICH THERE IS NO LONGER TURNOVER IN PRESENCE ABSENCE
# I.E., for replicate PCRs, how many counts can be observed of one OTU 
# where it is completely absent from another PCR?
turnover_thresholds <- no_turnover(otu_table[["scaled"]], metadata[["scaled"]][,colname_env_sample])
(turnover_threshold <- max(turnover_thresholds))
#-------------------------------------------------------------------------------

# THERE IS SOMETHING TO THIS:
# boxplot(apply(otu_table[["raw"]][1:4,1:20], 2, scale))

#-------------------------------------------------------------------------------
# Should OTUs be excluded that occur fewer than a threshold number of times?
EXCLUDE_RARE_OTUs <- TRUE
abundance_threshold <- turnover_threshold

# should the threshold be applied to each observation, or across all samples
# If at observation, any count < abundance threshold will be set to 0
# If across all samples, counts set to 0 for columns with sum < abundance_threshold

THRESHOLD_LEVEL <- "observation" # "observation" | "all_samples"

if(EXCLUDE_RARE_OTUs) {
  if(THRESHOLD_LEVEL == "all_samples"){
	  # how many OTUs will be retained?
	  sum(colSums(otu_table[["scaled"]]) > abundance_threshold)
	  # round(colSums(otu_scaled)/sum(otu_scaled)*100, digits = 2)
	  # plot(colSums(otu_scaled)/sum(otu_scaled), type = "l")
	  # exclude OTUs that were counted a total of fewer than 100 times across samples
	  otu_table[["filt"]] <- otu_table[["scaled"]][ ,
	      which(colSums(otu_table[["scaled"]]) > abundance_threshold)
	  ]
	  metadata[["filt"]] <- metadata[["scaled"]]
  } else if(THRESHOLD_LEVEL == "observation"){
  	# do it on a per-sample basis
    otu_table[["filt"]] <- otu_table[["scaled"]]
    otu_table[["filt"]][otu_table[["filt"]] < abundance_threshold] <- 0
    otu_table[["filt"]] <- strip_absent(otu_table[["filt"]])
    metadata[["filt"]]  <- metadata[["scaled"]]
  }
}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# take mean of OTU abundance across replicate PCRs
otu_table[["mean"]] <- do.call(rbind,
	lapply(
		split(as.data.frame(otu_table[["filt"]]), metadata[["filt"]][,colname_env_sample]),
	colMeans
	)
)
# reduce corresponding metdata
metadata[["mean"]] <- metadata[["filt"]][
   match(rownames(otu_table[["mean"]]), metadata[["filt"]][,colname_env_sample])
   , ]
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# take mean of unfiltered OTU abundance across replicate PCRs
otu_table[["mean_unfilt"]] <- do.call(rbind,
	lapply(
		split(as.data.frame(prop(otu_table[["clean"]])), metadata[["clean"]][,colname_env_sample]),
	colMeans
	)
)
# reduce corresponding metdata
metadata[["mean_unfilt"]] <- metadata[["clean"]][
   match(rownames(otu_table[["mean"]]), metadata[["clean"]][,colname_env_sample])
   , ]
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Log transform
otu_table[["log"]] <- log(otu_table[["mean"]] + 1)
metadata[["log"]]  <- metadata[["mean"]]
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# rescale each column to its maximum value, ranging from 0 to 1
otu_table[["scale01"]] <- apply(otu_table[["mean"]], MARGIN = 2, FUN = scale01)
metadata[["scale01"]]  <- metadata[["mean"]]
# par(mar = c(4,5,1,1))
# stripchart(as.data.frame(otu_table[["scale01"]][,1:20]), pch = 19, col = rgb(0,0,0, alpha = 0.2),
  # main = "", las = 1, xlab = "scaled proportion across samples")
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# exclude otus that have no spatial variance (occur in only one sample)
otu_table[["spvar"]] <- rm_single_rows(otu_table[["mean"]])
metadata[["spvar"]]  <- metadata[["mean"]]
#-------------------------------------------------------------------------------

if(!identical(length(otu_table), length(metadata))){
	warning("watch out, the number of OTU and metadata tables differs")
}
if(!identical(names(otu_table), names(metadata))){
	warning("watch out, the names of the OTU tables and metadata differ")
}
