# !!! SET WORKING DIRECTORY TO THIS PROJECT'S SUBDIRECTORY 'Analysis'
# Set directories from which to read/write data and write figures
analysis_dir <- dirname(file.choose())
setwd(analysis_dir)
data_dir <- file.path("..", "Data")
fig_dir <- file.path("..", "Figures")

# REQUIRES:
# 1. Working directory is set to the 'Analysis' folder within this project's folder
# 2. An OTU table (CSV) in the 'Data' folder
# 3. A metadata file (CSV) in the 'Data' folder containing information for all and only the samples in the OTU table

# LOAD OTU table
# path to otu file: rows are samples, columns are OTUs, cells are counts.
# rownames correspond to column "sample_id" in metadata
otu_table_filename <- "OTU_table.csv"

otu_table <- read.csv(
	file = file.path(data_dir, otu_table_filename), 
	row.names = 1, 
	stringsAsFactors = FALSE
	)
	
dim(otu_table)

# transpose the OTU table if "lib_" is found in the column names
# (this indicates a sample, which should be in rows)
if(length(grep("lib_", colnames(otu_table))) > 0){
	otu_table <- t(otu_table)
}
dim(otu_table)

# LOAD METADATA
# path to metadata file.
# column "sample_id" contains rownames of otu file.
metadata_filename <- "metadata_spatial.csv"

metadata <- read.csv(
	file = file.path(data_dir, metadata_filename), 
	stringsAsFactors = FALSE
	)

# exclude rows from OTU table that are not in metadata
otu_table <- otu_table[metadata$sample_id,]

# exclude OTUs not found in these samples
otu_table <- otu_table[,which(colSums(otu_table) > 0)]
dim(otu_table)
if( nrow(metadata) != nrow(otu_table) ) {
	warning("number of rows in metadata and otu table are not equal, but they should be.")
}

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


# exclude OTUs that were counted a total of fewer than 100 times across samples
otu_filt <- otu_scaled[,which(colSums(otu_scaled) > abundance_threshold)]
dim(otu_filt)

boxplot(otu_filt[,1:20])
# boxplot(log(otu_table[,1:20]))
boxplot(otu_table_prop[,1:20])
boxplot(scale(otu_table[,1:20]))

min(scale(otu_table[,1:20])+2)
boxplot(scale(otu_table[,1:20])+2)
