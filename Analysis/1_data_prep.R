
# holy macaroni, we need to decide on a standard orientation.
# "OTUs_fam_w30_top10.csv" is in rows=OTU, cols=SAMPLES
# YOU SHOULD PROBABLY DO THIS: start with script '0_load_data.R',
# use the resulting object 'otu_filt' in place of table_transposed
# OR use t(otu_filt) as object 'table_full'

# filename_orig <- file.path(data_dir, "OTUs_fam_w30_top10.csv")
# filename_out <- file.path(data_dir, "OTUs_top10.csv")

# table_full <- read.csv(filename_orig, row.names = 1)
# maybe ? table_full <- t(otu_filt)

# write the file
# write.csv(x = table_restricted, file = filename_out, quote = FALSE)

# TODO: append tables to a list like so:
# some_tables[["yet_another"]] <- matrix(data = rnorm(25), nrow = 5)

R_files <- list.files(path = "functions", pattern = "\\.R$", full.names = TRUE)
sapply(R_files, source)

otu_table_in <- otu_table_raw

# calculate proportional abundance of OTUs in each sample
otu_table_prop <- otu_table_in/rowSums(otu_table_in)


#-------------------------------------------------------------------------------
CHECK_FOR_OUTLIERS <- TRUE
# If you'd like to check for and remove inconsistent PCR replicates, go to:
if(CHECK_FOR_OUTLIERS) {
  cleaned <- find_bad_replicate(
  	my_table = otu_table_in,
  	my_metadata = metadata,
  	sample_id_column = colname_sampleid,
  	grouping_column = colname_env_sample, 
  	threshold_sd = 1, 
    save_pdf = TRUE, 
    pdf_path = file.path(fig_dir, "PCR_consistency.pdf"))
  otu_clean <- strip_absent(cleaned[[1]])
  metadata_clean <- cleaned[[2]]
  rm(cleaned)
}
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# RESCALE TO EQUAL SEQUENCING DEPTHS PER SAMPLE
RESCALE_SEQUENCING_DEPTH <- TRUE

if(RESCALE_SEQUENCING_DEPTH) {

	otu_scaled <- rescale_rowsums(otu_clean)

}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# CALCULATE THRESHOLD AT WHICH THERE IS NO LONGER TURNOVER IN PRESENCE ABSENCE
# I.E., for replicate PCRs, how many counts can be observed of one OTU 
# where it is completely absent from another PCR?
turnover_thresholds <- no_turnover(otu_scaled, metadata_clean[,colname_env_sample])
(turnover_threshold <- max(turnover_thresholds))
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# take mean of OTU abundance across replicate PCRs
otu_mean <- do.call(rbind,
	lapply(
		split(as.data.frame(otu_clean), metadata_clean[,colname_env_sample]),
	colMeans
	)
)
# reduce corresponding metdata
metadata_mean <- metadata_clean[match(rownames(otu_mean), metadata_clean[,colname_env_sample]),]
#-------------------------------------------------------------------------------

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
	  sum(colSums(otu_scaled) > abundance_threshold)
	  # round(colSums(otu_scaled)/sum(otu_scaled)*100, digits = 2)
	  # plot(colSums(otu_scaled)/sum(otu_scaled), type = "l")
	  # exclude OTUs that were counted a total of fewer than 100 times across samples
	  otu_filt <- otu_scaled[,which(colSums(otu_scaled) > abundance_threshold)]
	  dim(otu_filt)
  } else if(THRESHOLD_LEVEL == "observation"){
  	# do it on a per-sample basis
	otu_filt <- otu_mean
	otu_filt[otu_filt < abundance_threshold] <- 0
	otu_filt <- strip_absent(otu_filt)
	dim(otu_filt)
  }
}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Log transform
otu_log <- log(otu_filt + 1)
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# rescale each column to its maximum value, ranging from 0 to 1
otu_01 <- apply(otu_filt, MARGIN = 2, FUN = scale01)
par(mar = c(4,5,1,1))
stripchart(as.data.frame(otu_01[,1:20]), pch = 19, col = rgb(0,0,0, alpha = 0.2),
  main = "", las = 1, xlab = "scaled proportion across samples")
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# exclude otus that have no spatial variance (occur in only one sample)
otu_spvar <- rm_single_rows(otu_filt)
dim(otu_spvar)
#-------------------------------------------------------------------------------
