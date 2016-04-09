# Manipulate data in preparation for modeling in JAGS

data_dir <- file.path("..", "Data")
fig_dir <- file.path("..", "Figures")

# generate a file of only the samples in the metadata
# (assumes the metadata only has samples we're interested in rather than the full sequencing run)

# holy macaroni, we need to decide on a standard orientation.
# "OTUs_fam_w30_top10.csv" is in rows=OTU, cols=SAMPLES
# YOU SHOULD PROBABLY DO THIS: start with script '0_load_data.R',
# use the resulting object 'otu_filt' in place of table_transposed
# OR use t(otu_filt) as object 'table_full'

# filename_orig <- file.path(data_dir, "OTUs_fam_w30_top10.csv")

# filename_out <- file.path(data_dir, "OTUs_top10.csv")

# table_full <- read.csv(filename_orig, row.names = 1)
# maybe ? table_full <- t(otu_filt)

# load metadata
# metadata_file <- file.path(data_dir, "metadata_spatial.csv")
# metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)


# get stuff togther in a way we can use in STAN
sequenced_sample <- 1:nrow(metadata)
transect_position <- as.numeric(as.factor(metadata[,"dist_from_shore"]))
transect_line_rep <- as.numeric(as.factor(metadata[,"transect"]))
transect_line <- as.numeric(as.factor(paste(transect_line_rep, transect_position)))

# what is the name of the column containing the sample id (SEQUENCING samples)
colname_sampleid <- "sample_id"

metadata_exp <- cbind.data.frame(metadata, sequenced_sample, transect_position, transect_line)

# restrict to only samples in the metadata file
# !!! it will really screw things up if you don't convert the sample id column to a character vector!
table_restricted <- table_full[ , as.character(metadata[,colname_sampleid])]

# write the file
# write.csv(x = table_restricted, file = filename_out, quote = FALSE)

# COMBINE UP DATA
table_transposed <- t(table_restricted[,as.character(metadata[,colname_sampleid])])
# table_transposed <- otu_filt
if(
  identical(
    rownames(table_transposed),
    metadata[,colname_sampleid]
  )
){
  data_full <- cbind.data.frame(metadata_exp, table_transposed)
  if(identical(rownames(data_full), data_full[,colname_sampleid])){
    rownames(data_full) <- NULL
  } else {
    "something went wrong, check thyself before you wreck thyself"
  }
} else {
  warning("something does not match up, and you could screw up the data big time by proceeding.")
}

# data_full_path <- file.path(data_dir, "data_full_top10.csv")
data_full_path <- file.path(data_dir, "data_full_filt.csv")
# write.csv(x = data_full, file = data_full_path, row.names = FALSE)

# LONG FORMAT:
library(reshape2)
otu_long <- melt(otu_filt)
colnames(otu_long) <- c("sample_id", "OTU", "count")
otu_long$OTU <- as.character(otu_long$OTU)
otu_long$sample_id <- as.character(otu_long$sample_id)

row_match_long <- match(otu_long$sample_id , metadata_exp[, colname_sampleid])

data_full_long <- cbind.data.frame(metadata_exp[row_match_long,], otu_long[,c("OTU", "count")])

data_for_cast <- data.frame(
  sample_id = data_full_long$sample_id, 
  OTU = data_full_long$OTU, 
  count = data_full_long$count)

# note high value for sample "lib_B_tag_GCGCTC" in PCT-C-0500
# remove that sample?
remove_outlier <- TRUE
if(remove_outlier == TRUE){
  metadata <- metadata[metadata$sample_id != "lib_B_tag_GCGCTC",]
  metadata_exp <- metadata_exp[metadata_exp$sample_id != "lib_B_tag_GCGCTC",]
  otu_filt <- otu_filt[rownames(otu_filt) != "lib_B_tag_GCGCTC",]
  print("removed outlier PCR")
}


#-----------------------------------------------------------
# take mean of OTU abundance across replicate PCRs
otu_mean <- do.call(rbind, 
	lapply(
	lapply(
		split(otu_filt, metadata$env_sample_name), 
		matrix, ncol = ncol(otu_filt)
		), 
	colMeans
	)
)
colnames(otu_mean) <- colnames(otu_filt)


# I think I may have been crazy when I wrote the next few lines:
library(reshape2)
otu_long_mean <- melt(dcast(data = data_for_cast, formula = sample_id ~ OTU, fun.aggregate = mean, na.rm = TRUE))
colnames(otu_long_mean) <- c("sample_id", "OTU", "count")

row_match_mean <- match(otu_long_mean$sample_id , metadata_exp[, colname_sampleid])
otu_long_mean_full <- cbind.data.frame(metadata_exp[row_match_long,], otu_long[,c("OTU", "count")])

split(data_for_cast$count, f = c(data_for_cast$sample_id, data_for_cast$OTU))

# aggregate(x = data_for_cast, by = list(data_for_cast$OTU, data_for_cast$count), FUN = mean)
#-----------------------------------------------------------


data_full_long_path <- file.path(data_dir, "data_full_long.csv")
# write.csv(x = data_full_long, file = data_full_long_path, row.names = FALSE)




dcast()
dcast(data_full_long[1:400,], sample_id+OTU~count)
  data_full_long$sample_id

dcast(mtcars, )
head(mtcars)
# --------------- DATA TRIMMING ------------------------------------------
# RESTRICT TO ONLY A CERTAIN DISTANCE
# a table of counts of sequences (Z, length = ...)
counts_table_file <- file.path(data_dir, "OTUs_top10.csv")
counts_table <- read.csv(counts_table_file, row.names = 1)

# rows/rownames = samples, columns/colnames = taxa, cells = integer counts
# if table is incorrectly oriented, transpose it:
counts_table <- as.data.frame(t(counts_table))

# what is the name of the column containing distance from shore data 
colname_dist <- "dist_from_shore"

# what distance should we look at?
dist_of_interest <- 4000

# subset to only include samples taken at 4000m from shore
metadata <- droplevels(metadata[ metadata[ , colname_dist] == dist_of_interest , ])

# the ids of the sequenced samples we care about are given by
relevant_sampleids <- metadata[ , colname_sampleid]

# the counts that we care about are given by:
counts_table <- counts_table[relevant_sampleids , ]
write.csv(x = counts_table, file = "OTUs_top10_4000m.csv")
# ------------------------------------------------------------------------
