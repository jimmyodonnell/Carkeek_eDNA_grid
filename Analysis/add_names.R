# add organism names to an otu table
the_otu_table <- otu_mean

blast_file <- "query_hit_LCA.csv"

taxa_df <- read.csv(file.path(data_dir, blast_file))

head(taxa_df)

# remove abundance annotation from query sequence names
taxa_df$query_seq <- gsub(pattern = ";size=.*", replacement = "", taxa_df$query_seq)
head(taxa_df)


taxa_vector <- match(colnames(the_otu_table), taxa_df$query_seq)
nrow(taxa_df)

# which taxa are not in the blast file?
not_in_blastfile <- colnames(the_otu_table)[!colnames(the_otu_table) %in% taxa_df$query_seq]

# how abundant are they?
colSums(the_otu_table)[!colnames(the_otu_table) %in% taxa_df$query_seq]

# taxon ID of each otu
taxa_df$query_seq == 
colnames(the_otu_table)

colnames(taxa_df)
head(classification_df)

as.character(classification_df$query)

# columns in the taxa_df with taxon_id
as.character(taxa_df$LCA_id_all)
as.character(taxa_df$LCA_id_beste)

# column of the classification_df with taxon_id
as.character(classification_df$query)

match(as.character(taxa_df$LCA_id_all), as.character(classification_df$query))
match(as.character(classification_df$query), as.character(taxa_df$LCA_id_beste))

# what is the row of the taxa_df for each of the classification_df taxids
unique(match(as.character(classification_df$query), as.character(taxa_df$LCA_id_beste)))

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
  # out[i] <- (length(which(taxa_df$query_seq == seqnames[i])))
  print(taxa_df[which(taxa_df$query_seq == seqnames[i]), taxid_column] )
  taxids[i] <- taxa_df[which(taxa_df$query_seq == seqnames[i]), taxid_column]
  print(which(classification_df[,"query"] == taxids[i]))
  print(classification_list[[taxids[i]]])
  print(classification_list[[taxids[i]]][,])
}
as.character(classification_df[which(classification_df["rankname"] == "phylum"),"taxname"])
unique(as.character(classification_df[which(classification_df["rankname"] == "phylum"),"taxname"]))


classification_list[["6548"]]

family_beste <- taxa_df[taxa_vector, "family_beste"]
taxid_all <- taxa_df[taxa_vector, "LCA_id_all"]

# SET THE COLUMN OF THE TAXA TABLE HERE
names_alle <- gsub(" ", "_", as.character(taxa_df[taxa_vector, "LCA_name_all"] ))

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


gsub(" ", "_", as.character(taxa_df$LCA_name_beste[taxa_vector]))
otu_mean_LCT <- otu_mean

taxnames <- gsub(" ", "_", as.character(taxa_df$LCA_name_beste[taxa_vector]))
taxnames[is.na(taxnames)] <- "no_hit"

colnames(otu_mean_LCT) <- taxnames
