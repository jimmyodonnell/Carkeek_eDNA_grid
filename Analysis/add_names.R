# add organism names to an otu table
the_otu_table <- otu_filt

blast_file <- "query_hit_LCA.csv"

blast_df <- read.csv(
  file = file.path(data_dir, blast_file),
  stringsAsFactors = FALSE
)

head(blast_df)

# remove abundance annotation from query sequence names
blast_df$query_seq <- gsub(pattern = ";size=.*", replacement = "", blast_df$query_seq)
head(blast_df)


taxa_vector <- match(colnames(the_otu_table), blast_df$query_seq)
nrow(blast_df)

# which taxa are not in the blast file?
(not_in_blastfile <- colnames(the_otu_table)[
	!colnames(the_otu_table) %in% blast_df$query_seq
	])

# how abundant are they?
round(colSums(the_otu_table)[!colnames(the_otu_table) %in% blast_df$query_seq] / sum(the_otu_table), digits = 3)

# plot
# plot(
  # colSums(the_otu_table),
  # col = as.numeric(!colnames(the_otu_table) %in% blast_df$query_seq) + 1,
  # log = "y", ylab = "", las = 1
# )


taxon_column <- "family_all"
# LCA_name_all   | class_all   | order_all   | family_all
# LCA_name_beste | class_beste | order_beste | family_beste
table(blast_df$LCA_id_all   %in% classification_df$query)
table(blast_df$LCA_id_beste %in% classification_df$query)


table(blast_df$LCA_name_all   %in% classification_df$taxname)
table(blast_df$LCA_name_beste %in% classification_df$taxname)
table(is.na(blast_df$LCA_name_beste))
table(is.na(blast_df$LCA_name_beste))

# get the index where each of the blast results are found in the
# classification hierarchy data frame
blast_to_classif <- match(blast_df$LCA_name_beste, classification_df$taxname)

head(classification_df)
split(classification_df, classification_df$query)
classification_df$query[blast_to_classif]




taxa_by_otu <- blast_df[taxa_vector, taxon_column]
source("aggregate_cols.R")
taxon_table <- aggregate_cols(the_otu_table, col_groups = taxa_by_otu)


stop()
# OK pick up here

# columns in the blast_df with taxon_id
as.character(blast_df$LCA_id_all)
as.character(blast_df$LCA_id_beste)

# column of the classification_df with taxon_id
as.character(classification_df$query)

match(as.character(blast_df$LCA_id_all), as.character(classification_df$query))
match(as.character(classification_df$query), as.character(blast_df$LCA_id_beste))

# what is the row of the blast_df for each of the classification_df taxids
unique(match(as.character(classification_df$query), as.character(blast_df$LCA_id_beste)))

classification_list <- split(classification_df, classification_df$query)
classification_list[["1001339"]]
some_id <- "1001339"
classification_list[[some_id]]
taxid_column <- "LCA_id_beste"
seqnames <- colnames(otu_mean)
out <- vector()
taxids <- vector()
for(i in 1:length(seqnames)) {
  print(seqnames[i])
  # out[i] <- (length(which(blast_df$query_seq == seqnames[i])))
  print(blast_df[which(blast_df$query_seq == seqnames[i]), taxid_column] )
  taxids[i] <- blast_df[which(blast_df$query_seq == seqnames[i]), taxid_column]
  print(which(classification_df[,"query"] == taxids[i]))
  print(classification_list[[taxids[i]]])
  print(classification_list[[taxids[i]]][,])
}
as.character(classification_df[which(classification_df["rankname"] == "phylum"),"taxname"])
unique(as.character(classification_df[which(classification_df["rankname"] == "phylum"),"taxname"]))


classification_list[["6548"]]

family_beste <- blast_df[taxa_vector, "family_beste"]
taxid_all <- blast_df[taxa_vector, "LCA_id_all"]

# SET THE COLUMN OF THE TAXA TABLE HERE
names_alle <- gsub(" ", "_", as.character(blast_df[taxa_vector, "LCA_name_all"] ))

otu_naming <- otu_mean
colnames(otu_naming) <- names_alle

# using only column names (i.e. there should be duplicate column names):
aggregate_by_colname <- function(mat) {
  agg_mat <- do.call(cbind,
    lapply(unique(colnames(mat)), function(x) rowSums(as.matrix(mat[ , which(colnames(mat) == x)])))
  )
  colnames(agg_mat) <- unique(colnames(mat))
  return(agg_mat)
}

otu_named <- aggregate_by_colname(otu_naming)

if(sum(otu_named) != sum(otu_naming)){
	warning("watch it, you lost some OTUs when the data were collapsed by name")
}

# using the classifications data frame

head(classification_df)
classification_df

classification_list <- split(classification_df, classification_df$query)
classification_list[taxid_all[2]]


gsub(" ", "_", as.character(blast_df$LCA_name_beste[taxa_vector]))
otu_mean_LCT <- otu_mean

taxnames <- gsub(" ", "_", as.character(blast_df$LCA_name_beste[taxa_vector]))
taxnames[is.na(taxnames)] <- "no_hit"

colnames(otu_mean_LCT) <- taxnames
