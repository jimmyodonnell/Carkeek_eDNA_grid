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
