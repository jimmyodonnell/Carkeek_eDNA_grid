# !!! SET WORKING DIRECTORY TO THIS PROJECT'S SUBDIRECTORY 'Analysis'
# Set directories from which to read/write data and write figures
analysis_dir <- dirname(file.choose()) # choose this script file
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

# what is the name of the column containing the sample id (SEQUENCING samples)
colname_sampleid <- "sample_id"

# get stuff togther in a way we can use in STAN
sequenced_sample <- 1:nrow(metadata)
transect_position <- as.numeric(as.factor(metadata[,"dist_from_shore"]))
transect_line_rep <- as.numeric(as.factor(metadata[,"transect"]))
transect_line <- as.numeric(as.factor(paste(transect_line_rep, transect_position)))

metadata_exp <- cbind.data.frame(metadata, sequenced_sample, transect_position, transect_line)


# exclude rows from OTU table that are not in metadata
if(class(metadata[,colname_sampleid]) != "character"){
	warning("# !!! it will really screw things up if you don't convert the sample id column to a character vector!")
} else {
	otu_table <- otu_table[metadata[,colname_sampleid],]
}

# exclude OTUs not found in these samples
otu_table <- otu_table[,which(colSums(otu_table) > 0)]
dim(otu_table)
if( nrow(metadata) != nrow(otu_table) ) {
	warning("number of rows in metadata and otu table are not equal, but they should be.")
}

