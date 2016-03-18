data_dir <- file.path("..", "Data")
fig_dir <- file.path("..", "Figures")
# assemble metadata
metadata_file <- "metadata_spatial.csv"

coords_file <- "transects.csv"

coords_df <- read.csv(file.path(data_dir, coords_file), stringsAsFactors = FALSE)

sample_name <- paste("PCT", substr(coords_df$trans, 1, 1), sprintf("%04d", coords_df$point), sep = "-")
coords_df <- data.frame(coords_df, sample_name, stringsAsFactors = FALSE)[,c(1,2,3,4,5,7,6)]


metadata <- read.csv(file.path(data_dir, metadata_file), stringsAsFactors = FALSE)

coords_of_metadata <- match(metadata$env_sample_name, coords_df$sample_name)

metadata_with_coords <- cbind.data.frame(metadata, coords_df[coords_of_metadata,c("lon", "lat")], stringsAsFactors = FALSE)
write.csv(metadata_with_coords, file = "metadata_spatial.csv", quote = TRUE, row.names = FALSE)