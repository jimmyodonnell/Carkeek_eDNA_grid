
library(fpc)
library(vegan)


data_dir <- file.path("..", "Data")
fig_dir <- file.path("..", "Figures")


otu_table_filename <- "otu_table_filtered_per_samp.csv"
metadata_filename <- "metadata_spatial.csv"

otu_table <- read.csv(file.path(data_dir, otu_table_filename), row.names = 1)
metadata <- read.csv(file.path(data_dir, metadata_filename))


head(otu_table)[,1:10]

