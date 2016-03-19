# !!! SET WORKING DIRECTORY TO THIS PROJECT'S SUBDIRECTORY 'Analysis'
# Set directories from which to read/write data and write figures
data_dir <- file.path("..", "Data")
fig_dir <- file.path("..", "Figures")

# path to otu file: rows are samples, columns are OTUs, cells are counts.
# rownames correspond to column "sample_id" in metadata
otu_table_filename <- "otu_table_filtered_per_samp.csv"

# path to metadata file.
# column "sample_id" contains rownames of otu file.
metadata_filename <- "metadata_spatial.csv"

otu_table <- read.csv(file.path(data_dir, otu_table_filename), row.names = 1, stringsAsFactors = FALSE)
metadata <- read.csv(file.path(data_dir, metadata_filename), stringsAsFactors = FALSE)
dist_along_shore <- c(1000,2000,0)[as.numeric(as.factor(metadata$transect))]
metadata <- cbind.data.frame(metadata, dist_along_shore, stringsAsFactors = FALSE)

# exclude rows from OTU table that are not in metadata
otu_table <- otu_table[metadata$sample_id,]

otu_table_prop <- otu_table/rowSums(otu_table)

head(otu_table)[,1:10]
head(metadata)
metadata$sample_id


boxplot(otu_table[,1:20])
# boxplot(log(otu_table[,1:20]))
boxplot(otu_table_prop[,1:20])
boxplot(scale(otu_table[,1:20]))

min(scale(otu_table[,1:20])+2)
boxplot(scale(otu_table[,1:20])+2)
