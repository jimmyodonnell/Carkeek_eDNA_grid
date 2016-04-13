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

# this function rescales a numeric vector to 0 and 1
scale01 <- function(x){(x-min(x))/(max(x)-min(x))}

# RESCALE TO EQUAL SEQUENCING DEPTHS PER SAMPLE
# calculate the minimum number of reads assigned to these OTUs in these samples
minreads <- min(rowSums(otu_table))

# calculate proportional abundance of OTUs in each sample
otu_table_prop <- otu_table/rowSums(otu_table)

# scale the proportional abundance of the OTU in each sample to the minimum number of reads
otu_scaled <- otu_table_prop * minreads

# ignore counts OTUs found fewer than 0.5 times (anything greater would be rounded to 1)
otu_scaled[otu_scaled < 0.5] <- 0

# round to whole numbers
otu_scaled <- round(otu_scaled)
dim(otu_scaled)

# again exclude OTUs not found in these samples
otu_scaled <- otu_scaled[,which(colSums(otu_scaled) > 0)]
dim(otu_scaled)

# to re-order by the abundance in THESE samples (i.e. not those from samples from elsewhere)
otu_scaled <- otu_scaled[,order(colSums(otu_scaled), decreasing = TRUE)]
otu_scaled[1:10,1:10]

# how many OTUs occur more than a threshold number of times?
abundance_threshold <- 100
sum(colSums(otu_scaled) > abundance_threshold)
# round(colSums(otu_scaled)/sum(otu_scaled)*100, digits = 2)
# plot(colSums(otu_scaled)/sum(otu_scaled), type = "l")

# exclude OTUs that were counted a total of fewer than 100 times across samples
otu_filt <- otu_scaled[,which(colSums(otu_scaled) > abundance_threshold)]
dim(otu_filt)

boxplot(otu_filt[,1:20])
# boxplot(log(otu_table[,1:20]))
boxplot(otu_table_prop[,1:20])
boxplot(scale(otu_table[,1:20]))


#-------------------------------------------------------------------------------
# CHECK FOR OUTLIERS
# If you'd like to check for and remove replicates that seem inconsistent, go to:
# 'dissimilarity.R'
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# take mean of OTU abundance across replicate PCRs
otu_mean <- do.call(rbind, 
	lapply(
		split(as.data.frame(otu_filt), metadata[,colname_env_sample]), 
	colMeans
	)
)
metadata_mean <- metadata[match(rownames(otu_mean), metadata[,colname_env_sample]),]

# I think I may have been crazy when I wrote the next few lines:
library(reshape2)
otu_long_mean <- melt(dcast(data = data_for_cast, formula = sample_id ~ OTU, fun.aggregate = mean, na.rm = TRUE))
colnames(otu_long_mean) <- c("sample_id", "OTU", "count")
row_match_mean <- match(otu_long_mean$sample_id , metadata_exp[, colname_sampleid])
otu_long_mean_full <- cbind.data.frame(metadata_exp[row_match_long,], otu_long[,c("OTU", "count")])
split(data_for_cast$count, f = c(data_for_cast$sample_id, data_for_cast$OTU))
# aggregate(x = data_for_cast, by = list(data_for_cast$OTU, data_for_cast$count), FUN = mean)
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# rescale each column to its maximum value, ranging from 0 to 1
otu_01 <- apply(otu_mean, MARGIN = 2, FUN = scale01)
stripchart(as.data.frame(otu_01[,1:20]), pch = 19, col = rgb(0,0,0, alpha = 0.2), las = 1)
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# COMBINE UP DATA
if(
  identical(
    rownames(otu_table),
    metadata[,colname_sampleid]
  )
){
  data_full <- cbind.data.frame(metadata_exp, otu_table)
  if(identical(rownames(data_full), data_full[,colname_sampleid])){
    rownames(data_full) <- NULL
  } else {
    "something went wrong, check thyself before you wreck thyself"
  }
} else {
  warning("something does not match up, and you could screw up the data big time by proceeding.")
}
#-------------------------------------------------------------------------------

# data_full_path <- file.path(data_dir, "data_full_top10.csv")
data_full_path <- file.path(data_dir, "data_full_filt.csv")
# write.csv(x = data_full, file = data_full_path, row.names = FALSE)

#-------------------------------------------------------------------------------
# LONG FORMAT:
library(reshape2) # melt()
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




data_full_long_path <- file.path(data_dir, "data_full_long.csv")
# write.csv(x = data_full_long, file = data_full_long_path, row.names = FALSE)

