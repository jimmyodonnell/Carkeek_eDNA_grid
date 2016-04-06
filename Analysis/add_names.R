# add organism names

blast_file <- "query_hit_LCA.csv"

taxa_df <- read.csv(file.path(data_dir, blast_file))

head(taxa_df)

# remove abundance annotation from query sequence names

taxa_df$query_seq <- gsub(pattern = ";size=.*", replacement = "", taxa_df$query_seq)
head(taxa_df)


taxa_vector <- match(colnames(otu_filt), taxa_df$query_seq)

as.character(taxa_df$class_beste[taxa_vector])
