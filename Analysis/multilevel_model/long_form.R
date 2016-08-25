# Prep data and manipulate for multilevel modeling
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
