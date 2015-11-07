# put OTU data in long format:

library(reshape2)

setwd("/Users/jimmy.odonnell/Projects/Carkeek_eDNA_grid/Analysis")

data_dir <- file.path("..", "Data")
fig_dir <- file.path("..", "Figures")


# a table of counts of sequences (Z, length = ...)
counts_file_path <- file.path(data_dir, "OTUs_top10_4000m.csv")

counts_table <- read.csv(file = counts_file_path, row.names = 1)
counts_table <- counts_table[,1:10] # use only the top ten for now

# rows/rownames = samples, columns/colnames = taxa, cells = integer counts
# if table is incorrectly oriented, transpose it:
# counts_table <- as.data.frame(t(counts_table))

# REMEMBER TO REMOVE CONTROL TAXON IF IT'S STILL THERE!
# counts_table <- counts_table[,! names(counts_table) %in% "DUP_3"]


# load metadata
metadata_file_path <- file.path(data_dir, "metadata_spatial.csv")
metadata <- read.csv(file = metadata_file_path, stringsAsFactors = FALSE)

# what is the name of the column containing the sample id (SEQUENCING samples)
colname_sampleid <- "sample_id"

# what is the name of the column containing PCR replicate levels
colname_pcr <- "PCR_replicate"

# what is the name of the column containing position in x direction (along shore) replicate levels
colname_posx <- "transect"

# what is the name of the column containing position in y direction (from shore) replicate levels
colname_posy <- "dist_from_shore"

# give a name to call the variable describing the things that were counted
varname_taxa <- "taxon"

# give a name to call the variable describing what the counts are
varname_counts <- "reads"


# of each taxon (length = I)
taxa <- colnames(counts_table)

# from each PCR replicate (length = J)
pcr <- as.character(unique(metadata[, colname_pcr]))

# at each position in x direction (length = K)
posx <- unique(metadata[, colname_posx])

# at each position in y direction (length = L)
posy <- as.character(unique(metadata[, colname_posy]))

# *** NOTE *** #
# for the automated array approach to work, rows of metadata and otu table must line up!
if(nrow(counts_table) != nrow(metadata)){
	if( nrow(counts_table) > nrow(metadata) ){
		# extract and order the relevant rows of the OTU table
		counts_rel <- counts_table[metadata[, colname_sampleid], ]
		metadata_rel <- metadata
		# counts_table <- counts_rel
	} else if( nrow(counts_table) < nrow(metadata) ){
		# extract only the relevant rows of metadata
		metadata_rel <- metadata[match(rownames(counts_table), metadata[,colname_sampleid]),]
		counts_rel <- counts_table
		# metadata <- metadata_rel
	}
} else if(nrow(counts_table) == nrow(metadata)){
	counts_rel <- counts_table
	metadata_rel <- metadata
}

# match(metadata[, colname_sampleid], rownames(counts_table))

# make sure the order of the samples in the metadata and OTU table are the same
identical(
	rownames(counts_rel), 
	metadata_rel[,colname_sampleid]
	)

# reshape into a long-format dataframe
data_l <- reshape(cbind(metadata_rel, counts_rel), 
  varying = colnames(counts_rel), # aka taxa
  v.names = varname_counts,
  timevar = varname_taxa, 
  times = colnames(counts_rel), 
  new.row.names = 1:(nrow(counts_rel)*ncol(counts_rel)),
  direction = "long")[c(colname_sampleid, colname_pcr, colname_posx, colname_posy, varname_taxa, varname_counts)]



################################################################################
# THREE DIMENSIONAL ARRAY:

# this reorders the names of taxa
data_array <- acast(
				data = data_l, 
				formula =  list(colname_pcr, varname_taxa, colname_posx), 
				value.var = varname_counts
				)

# ... and thus this fucks up the names 
# dimnames(data_array) <- list(
	# pcr = pcr,
	# taxa = taxa,
	# posx = posx
	# )

# so that if you're comparing this to the original way I created the array, this fails:
identical(counts, data_array)

# but this fixes it:
data_array <- data_array[,taxa,]
dimnames(data_array) <- list(
	pcr = pcr,
	taxa = taxa,
	posx = posx
	)
identical(counts, data_array)

data_array[ pcr[4] ,         ,         ]
data_array[        , taxa[8] ,         ]
data_array[        ,         , posx[2] ]



################################################################################
# FOUR DIMENSIONAL ARRAY:

counts_file_path <- file.path(data_dir, "otu_table_filtered_per_samp.csv")

# this reorders the names of taxa
data_array <- acast(
				data = data_l, 
				formula =  list(colname_pcr, varname_taxa, colname_posx, colname_posy), 
				value.var = varname_counts
				)

# but this fixes it:
data_array <- data_array[,taxa,,]
dimnames(data_array) <- list(
	pcr = pcr,
	taxa = taxa,
	posx = posx, 
	posy = posy
	)
identical(counts, data_array)

data_array[ pcr[4] ,         ,         ,         ]
data_array[        , taxa[9] ,         ,         ]
data_array[        ,         , posx[3] ,         ]
data_array[        ,         ,         , posy[1] ]


