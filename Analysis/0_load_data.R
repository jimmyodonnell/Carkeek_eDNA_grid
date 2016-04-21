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

#-------------------------------------------------------------------------------
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

# what is the name of the column containing the sample id from which subsamples were drawn (environmental samples)
colname_env_sample <- "env_sample_name"

# what is the name of the column containing the adjusted x and y coordinates
colname_xcoord <- "dist_along_shore"
colname_ycoord <- "dist_from_shore"

colname_lon <- "lon"
colname_lat <- "lat"

# add columns that will be compatible with
sequenced_sample <- 1:nrow(metadata)
transect_position <- as.numeric(as.factor(metadata[,colname_ycoord]))
transect_line_rep <- as.numeric(as.factor(metadata[,colname_xcoord])) #metadata[,"transect"]
transect_line <- as.numeric(as.factor(paste(transect_line_rep, transect_position)))

metadata <- cbind.data.frame(metadata, sequenced_sample, transect_position, transect_line)


#-------------------------------------------------------------------------------
# LOAD OTU table
# path to otu file: rows are samples, columns are OTUs, cells are counts.
# rownames correspond to column "sample_id" in metadata
otu_table_filename <- "OTU_table.csv"

otu_table_raw <- read.csv(
	file = file.path(data_dir, otu_table_filename),
	row.names = 1,
	stringsAsFactors = FALSE
	)

USE_SODM_TABLE <- FALSE
if (USE_SODM_TABLE) {
	otu_sodm_file <- "SODM/OTUs_BayesianVetted_OTUs.csv"
	otu_table_raw <- read.csv(
		file = otu_sodm_file,
		row.names = 1,
		stringsAsFactors = FALSE
		)

}
dim(otu_table_raw)

#-------------------------------------------------------------------------------
# transpose the OTU table if "lib_" is found in the column names
# (this indicates a sample, which should be in rows)
if(any(metadata[,colname_sampleid] %in% colnames(otu_table_raw))){
	otu_table_raw <- t(otu_table_raw)
}
dim(otu_table_raw)

#-------------------------------------------------------------------------------
# exclude rows from OTU table that are not in metadata
if(class(metadata[,colname_sampleid]) != "character"){
	warning("# !!! it will really screw things up if you don't convert the sample id column to a character vector!")
} else {
	otu_table_raw <- otu_table_raw[metadata[,colname_sampleid],]
}

#-------------------------------------------------------------------------------
# exclude OTUs not found in these samples
otu_table_raw <- otu_table_raw[,which(colSums(otu_table_raw) > 0)]
dim(otu_table_raw)
if( nrow(metadata) != nrow(otu_table_raw) ) {
	warning("number of rows in metadata and otu table are not equal, but they should be.")
}
