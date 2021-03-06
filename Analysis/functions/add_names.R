add_names <- function(
# add organism names to an otu table
# requires:
  otu_table = otu_table[["mean"]], 
  blast_file = file.path(data_dir, "query_hit_LCA.csv"), 
  taxon_column = "LCA_name_beste"
){

# specify which column to get taxon names from
# other taxon columns to use:
# LCA_name_all   | class_all   | order_all   | family_all
# LCA_name_beste | class_beste | order_beste | family_beste

the_table <- otu_table

blast_df <- read.csv(
  file = blast_file,
  stringsAsFactors = FALSE
)
head(blast_df)

# remove abundance annotation from query sequence names
blast_df$query_seq <- gsub(
  pattern = ";size=.*", replacement = "", blast_df$query_seq)
# head(blast_df)
# nrow(blast_df)

# make a vector to link OTUs to their corresponding row in blast_df
taxa_vector <- match(colnames(the_table), blast_df$query_seq)

# which taxa are not in the blast file?
not_in_blastfile <- colnames(the_table)[
	!colnames(the_table) %in% blast_df$query_seq
	]

# how abundant are they, in percent of the total reads
percents <- round(colSums(the_table)[!colnames(the_table) %in% blast_df$query_seq] / sum(the_table), digits = 3)*100
percent_unaccounted <- data.frame(
  otu_name   = names(percents),
  percentage = paste(as.character(percents), "%", sep = "")
)
warning_msg <- "The following OTUs are missing from the blast results:"

for(i in list(warning_msg, percent_unaccounted)){print(i)}

# plot
# plot(
  # colSums(the_table),
  # col = as.numeric(!colnames(the_table) %in% blast_df$query_seq) + 1,
  # log = "y", ylab = "", las = 1
# )

taxa_by_otu <- blast_df[taxa_vector, taxon_column]
taxon_table <- aggregate_cols(the_table, col_groups = taxa_by_otu)
return(taxon_table)

}

# GO TO life_history.R


# columns in the blast_df with taxon_id
# as.character(blast_df$LCA_id_all)
# as.character(blast_df$LCA_id_beste)

# column of the classification_df with taxon_id
# as.character(classification_df$query)

# match(as.character(blast_df$LCA_id_all), as.character(classification_df$query))
# match(as.character(classification_df$query), as.character(blast_df$LCA_id_beste))

# what is the row of the blast_df for each of the classification_df taxids
# unique(match(as.character(classification_df$query), as.character(blast_df$LCA_id_beste)))

# family_beste <- blast_df[taxa_vector, "family_beste"]
# taxid_all <- blast_df[taxa_vector, "LCA_id_all"]

# # SET THE COLUMN OF THE TAXA TABLE HERE
# names_alle <- gsub(" ", "_", as.character(blast_df[taxa_vector, "LCA_name_all"] ))
# colnames(otu_naming) <- names_alle


# if(sum(otu_named) != sum(otu_naming)){
	# warning("watch it, you lost some OTUs when the data were collapsed by name")
# }

# # using the classifications data frame

# gsub(" ", "_", as.character(blast_df$LCA_name_beste[taxa_vector]))
# otu_mean_LCT <- otu_mean

# taxnames <- gsub(" ", "_", as.character(blast_df$LCA_name_beste[taxa_vector]))
# taxnames[is.na(taxnames)] <- "no_hit"

# colnames(otu_mean_LCT) <- taxnames
