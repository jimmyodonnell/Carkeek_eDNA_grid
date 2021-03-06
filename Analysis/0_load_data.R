# !!! SET WORKING DIRECTORY TO THIS PROJECT'S SUBDIRECTORY 'Analysis'
# Set directories from which to read/write data and write figures
INTERACTIVE <- FALSE
if(INTERACTIVE){
  analysis_dir <- dirname(file.choose()) # choose this script file
  setwd(analysis_dir)
} else {
  analysis_dir <- getwd()
}
data_dir <- file.path("..", "Data")
fig_dir <- file.path("..", "Figures")

# REQUIRES:
# 1. Working directory is set to the 'Analysis' folder within this project's folder
# 2. An OTU table (CSV) in the 'Data' folder
# 3. A metadata file (CSV) in the 'Data' folder containing information for all and only the samples in the OTU table


#-------------------------------------------------------------------------------
# LOAD FUNCTIONS
R_files <- list.files(path = "functions", pattern = "\\.R$", full.names = TRUE)
sapply(R_files, source)
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# LOAD PACKAGES
required_packages <- c(
  'ape', 
  'cluster', 
  'colorspace', 
  'dplyr', 
  'fpc', 
  'geoR', 
  'geosphere', 
  'ggplot2', 
  'mgcv', 
  'propagate', 
  'raster', 
  'rasterVis', 
  'reshape2', 
  'vegan'
)
sapply(required_packages, load_package)
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# LOAD METADATA
# path to metadata file.
# column "sample_id" contains rownames of otu file.
metadata_filename <- "metadata_spatial.csv"

metadata <- list()

metadata[["raw"]] <- read.csv(
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
sequenced_sample <- 1:nrow(metadata[["raw"]])
transect_position <- as.numeric(as.factor(metadata[["raw"]][,colname_ycoord]))
transect_line_rep <- as.numeric(as.factor(metadata[["raw"]][,colname_xcoord])) #metadata[,"transect"]
transect_line <- as.numeric(as.factor(paste(transect_line_rep, transect_position)))

metadata[["raw"]] <- cbind.data.frame(
    metadata[["raw"]], sequenced_sample, transect_position, transect_line)


#-------------------------------------------------------------------------------
# LOAD OTU table
# path to otu file: rows are samples, columns are OTUs, cells are counts.
# rownames correspond to column "sample_id" in metadata
otu_table_filename <- "OTU_table.csv"

otu_table          <- list()

USE_SODM_TABLE <- FALSE
if (USE_SODM_TABLE) {
	otu_sodm_file <- "SODM/OTUs_BayesianVetted_OTUs.csv"
	otu_table[["raw"]] <- read.csv(
		file = otu_sodm_file,
		row.names = 1,
		stringsAsFactors = FALSE
	)
} else {
	otu_table[["raw"]] <- read.csv(
		file = file.path(data_dir, otu_table_filename),
		row.names = 1,
		stringsAsFactors = FALSE
	)
}
dim(otu_table[["raw"]])

#-------------------------------------------------------------------------------
# transpose the OTU table if "lib_" is found in the column names
# (this indicates a sample, which should be in rows)
if(any(metadata[["raw"]][,colname_sampleid] %in% colnames(otu_table[["raw"]]))){
	otu_table[["raw"]] <- t(otu_table[["raw"]])
}
dim(otu_table[["raw"]])

#-------------------------------------------------------------------------------
# exclude rows from OTU table that are not in metadata and correct order
if(class(metadata[["raw"]][,colname_sampleid]) != "character"){
	warning("# !!! it will really screw things up if you don't convert the sample id column to a character vector!")
} else {
	otu_table[["raw"]] <- otu_table[["raw"]][metadata[["raw"]][,colname_sampleid],]
}

#-------------------------------------------------------------------------------
# exclude OTUs not found in any of these samples
otu_table[["raw"]] <- otu_table[["raw"]][,which(colSums(otu_table[["raw"]]) > 0)]
dim(otu_table[["raw"]])
if( nrow(metadata[["raw"]]) != nrow(otu_table[["raw"]]) ) {
	warning("number of rows in metadata and otu table are not equal, but they should be.")
}
