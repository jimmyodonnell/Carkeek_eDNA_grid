# Prep data and manipulate for multilevel modeling

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

# this function rescales a numeric vector to 0 and 1
scale01 <- function (x) {(x-min(x))/(max(x)-min(x))}

# this function rescales a matrix to contain only 0 and 1
as.binary <- function (a_matrix) {
	bin_mat <- a_matrix
	bin_mat[bin_mat > 0] <- 1
	return(bin_mat)
}

present_in_all_rows <- function(x)
{
# this function removes columns which are not > 0 in all rows
	return(x[ , colSums(x > 0) >= nrow(x)])
}

strip_absent <- function(x)
{
# this function removes any columns for which the sum is <= 0
	return(x[,which(colSums(x) > 0)])
}

otu_table_in <- otu_table_raw

#-------------------------------------------------------------------------------
# RESCALE TO EQUAL SEQUENCING DEPTHS PER SAMPLE
RESCALE_SEQUENCING_DEPTH <- TRUE

if(RESCALE_SEQUENCING_DEPTH) {

	# calculate the minimum number of reads assigned to these OTUs in these samples
	minreads <- min(rowSums(otu_table_in))

	# calculate proportional abundance of OTUs in each sample
	otu_table_prop <- otu_table_in/rowSums(otu_table_in)

	# scale the proportional abundance of the OTU in each sample to the minimum number of reads
	otu_scaled <- otu_table_prop * minreads

	# ignore counts OTUs found fewer than 0.5 times (anything greater would be rounded to 1)
	otu_scaled[otu_scaled < 0.5] <- 0

	# round to whole numbers
	otu_scaled <- round(otu_scaled)
	dim(otu_scaled)

	# again exclude OTUs not found in these samples
	otu_scaled <- strip_absent(otu_scaled)
	dim(otu_scaled)

	# to re-order by the abundance in THESE samples (i.e. not those from samples from elsewhere)
	otu_scaled <- otu_scaled[,order(colSums(otu_scaled), decreasing = TRUE)]
	otu_scaled[1:10,1:10]

}


#-------------------------------------------------------------------------------
# Should OTUs be excluded that occur fewer than a threshold number of times?
EXCLUDE_RARE_OTUs <- TRUE
abundance_threshold <- 4

if(EXCLUDE_RARE_OTUs) {
  # how many OTUs will be retained?
  sum(colSums(otu_scaled) > abundance_threshold)
  # round(colSums(otu_scaled)/sum(otu_scaled)*100, digits = 2)
  # plot(colSums(otu_scaled)/sum(otu_scaled), type = "l")
  # exclude OTUs that were counted a total of fewer than 100 times across samples
  otu_filt <- otu_scaled[,which(colSums(otu_scaled) > abundance_threshold)]
  dim(otu_filt)
}


# boxplot(otu_filt[,1:20])
# boxplot(log(otu_table_in[,1:20]))
# boxplot(otu_table_prop[,1:20])
# boxplot(scale(otu_table_in[,1:20]))


#-------------------------------------------------------------------------------
CHECK_FOR_OUTLIERS <- TRUE
# If you'd like to check for and remove inconsistent PCR replicates, go to:
if(CHECK_FOR_OUTLIERS) {
  source('dissimilarity.R')
  cleaned <- find_bad_PCR(
  	my_table = otu_filt,
  	my_metadata = metadata,
  	sample_id_column = colname_sampleid,
  	grouping_column = colname_env_sample)
  otu_clean <- cleaned[[1]]
  metadata_clean <- cleaned[[2]]
  rm(cleaned)
}
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
# rescale each column to its maximum value, ranging from 0 to 1
otu_01 <- apply(otu_mean, MARGIN = 2, FUN = scale01)
par(mar = c(4,5,1,1))
stripchart(as.data.frame(otu_01[,1:20]), pch = 19, col = rgb(0,0,0, alpha = 0.2),
  main = "", las = 1)
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# exclude otus that have no spatial variance (occur in only one sample)
source("rm_single_rows.R")
otu_spvar <- rm_single_rows(otu_mean)
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# COMBINE UP DATA
if(
  identical(
    rownames(otu_table_in),
    metadata[,colname_sampleid]
  )
){
  data_full <- cbind.data.frame(metadata, otu_table_in)
  if(identical(rownames(data_full), data_full[,colname_sampleid])){
    rownames(data_full) <- NULL
  } else {
    "something went wrong, check thyself before you wreck thyself"
  }
} else {
  warning("something does not match up, and you could screw up the data big time by proceeding.")
}

# data_full_path <- file.path(data_dir, "data_full_filt.csv") # "data_full_top10.csv"
# write.csv(x = data_full, file = data_full_path, row.names = FALSE)
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# LONG FORMAT:
library(reshape2) # melt()
otu_long <- melt(otu_clean)
colnames(otu_long) <- c("sample_id", "OTU", "count")
otu_long$OTU <- as.character(otu_long$OTU)
otu_long$sample_id <- as.character(otu_long$sample_id)

row_match_long <- match(otu_long$sample_id , metadata[, colname_sampleid])

data_full_long <- cbind.data.frame(metadata[row_match_long,], otu_long[,c("OTU", "count")])

data_full_long_path <- file.path(data_dir, "data_full_long.csv")
# write.csv(x = data_full_long, file = data_full_long_path, row.names = FALSE)
# I think I may have been crazy when I wrote the next few lines:
library(reshape2)
data_for_cast <- data.frame(
  sample_id = data_full_long$sample_id,
  OTU = data_full_long$OTU,
  count = data_full_long$count)

otu_long_mean <- melt(dcast(data = data_for_cast, formula = sample_id ~ OTU, fun.aggregate = mean, na.rm = TRUE))
colnames(otu_long_mean) <- c("sample_id", "OTU", "count")
row_match_mean <- match(otu_long_mean$sample_id , metadata[, colname_sampleid])
otu_long_mean_full <- cbind.data.frame(metadata[row_match_long,], otu_long[,c("OTU", "count")])
split(data_for_cast$count, f = c(data_for_cast$sample_id, data_for_cast$OTU))
# aggregate(x = data_for_cast, by = list(data_for_cast$OTU, data_for_cast$count), FUN = mean)
