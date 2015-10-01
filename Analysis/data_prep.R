# Manipulate data in preparation for modeling in JAGS

filename_orig <- "/Users/threeprime/Documents/GoogleDrive/Kelly_Lab/Projects/Carkeek_eDNA_grid/Data/OTUs_fam_w30_top10.csv"

filename_out <- file.path(dirname(filename_orig), "top10_OTUs.csv")

metadata <- read.csv("/Users/threeprime/Documents/GoogleDrive/Kelly_Lab/Projects/Carkeek_eDNA_grid/Data/metadata_spatial.csv")

column_name <- "sample_id"

table_orig <- read.csv(filename_orig, row.names = 1)


# restrict to only samples in the metadata file
# it will really screw things up if you don't convert the sample id column to a character vector!
# table_new <- table_orig[,metadata$sample_id]
table_new <- table_orig[,as.character(metadata$sample_id)]


# write the file
write.csv(x = table_new, file = filename_out, quote = FALSE)


# --------------- DATA TRIMMING ------------------------------------------
# a table of counts of sequences (Z, length = ...)
counts_table <- read.csv("/Users/threeprime/Documents/GoogleDrive/Kelly_Lab/Projects/Carkeek_eDNA_grid/Data/OTUs_top10.csv", row.names = 1)
# rows/rownames = samples, columns/colnames = taxa, cells = integer counts
# if table is incorrectly oriented, transpose it:
counts_table <- as.data.frame(t(counts_table))

# load metadata
metadata <- read.csv("/Users/threeprime/Documents/GoogleDrive/Kelly_Lab/Projects/Carkeek_eDNA_grid/Data/metadata_spatial.csv", stringsAsFactors = FALSE)

# what is the name of the column containing the sample id (SEQUENCING samples)
colname_sampleid <- "sample_id"

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
