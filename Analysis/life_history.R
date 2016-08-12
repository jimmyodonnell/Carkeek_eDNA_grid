# add life history data

life_history_file <- "life_history.csv"
life_history <- read.csv(
  file = file.path(data_dir, life_history_file),
  stringsAsFactors = FALSE
)

head(life_history)
tail(life_history)

# Exclude data for groups where no data is given but a higher level is
life_history <- life_history[
!apply(life_history, 1, function(x) all(is.na(x[5:8])))
,]
# lh_taxon <- "taxon_name"

# join sequence count table and life history table, using blast results and classifications as intermediate

classification_file <- "classification_2016-04-19T125449.csv"

classification_df <- read.csv(
	file = file.path(data_dir, classification_file),
	row.names = 1,
	stringsAsFactors = FALSE
)
head(classification_df)
classification_col <- "taxname"


# which of the taxa in the classification df have data?
# classification_df$taxname %in% life_history$taxon_name

# Add a column to the classification_df, specifiying the corresponding 
# row in the life history data
classification_df$lifehist_row <- match(
	classification_df$taxname, life_history$taxon_name
)
head(classification_df)


##########
# How many of the values from a given column of the blast_df are in classification_df
table(blast_df$LCA_id_all   %in% classification_df$query)
table(blast_df$LCA_id_beste %in% classification_df$query)
table(blast_df$LCA_name_all   %in% classification_df$taxname) # all 1153
table(blast_df$LCA_name_beste %in% classification_df$taxname) # all 1153

# table(is.na(blast_df$LCA_name_beste))
# table(is.na(blast_df$LCA_name_beste))

# get the index where each of the blast results are found in the
# classification hierarchy data frame
blast_to_classif <- match(blast_df$LCA_name_beste, classification_df$taxname)

# I don't think I need to do the split anymore.
# classification_by_query <- split(classification_df, classification_df$query)
# has_numbers <- function(x){
	# any(!is.na(x[,"lifehist_row"]))
# }
# table(sapply(classification_by_query, has_numbers))

# get the first occurrence of a taxon in the classification_df
get_lh_data <- function(taxon){
	if(class(taxon) != "character"){
		stop("input should be a character string specifying a taxon name")
	}
	# show all of the rows of the query containing the taxon
	temp <- classification_df[
		classification_df$query == classification_df$query[
			match(taxon, classification_df$taxname)],
	]
	has_data <- temp$lifehist_row[1:match(taxon, temp$taxname)]
	return(temp$lifehist_row[max(which(!is.na(has_data)))])
}
get_lh_data(1)
get_lh_data("Hominidae")

blast_df$lifehist_row <- sapply(blast_df$LCA_name_beste, get_lh_data)

life_history[blast_df$lifehist_row,]
otu_to_lifehist <- blast_df[taxa_vector, "lifehist_row"]
life_history[otu_to_lifehist,"range_gamete_km"]
head(life_history)

life_history_OTU <- data.frame(
  query_taxon = taxa_by_otu, 
  life_history[otu_to_lifehist,], 
  row.names = NULL
)

not_marine <- which(
# exclude terrestrial
  life_history_OTU$adult_habitat == "terrestrial" |
# exclude freshwater
  life_history_OTU$adult_habitat == "freshwater"
)
the_otu_table[,-not_marine]
life_history_OTU[is.na(not_marine), ]
# not_marine[]

# gets rows meeting certain conditions:
low_dispersal <- 	life_history_OTU$adult_habitat != "terrestrial" &
	life_history_OTU$adult_habitat != "freshwater" &
	life_history_OTU$range_gamete_km < 10 &
	life_history_OTU$range_larva_km < 10 &
	life_history_OTU$range_adult_km < 10
low_dispersal[is.na(low_dispersal)] <- FALSE
life_history_OTU[ low_dispersal, ]
the_otu_table[,which(low_dispersal)]

# what percent of sequences belonging to organisms with high dispersal potential
round(sum(the_otu_table[,which(!low_dispersal)])/sum(the_otu_table) * 100, 1)

# what percent of otus were matched to organisms with high dispersal potential
round(sum(!low_dispersal)/length(low_dispersal) * 100, 1)

otu_sum_by_gamete_range <- split(colSums(the_otu_table), life_history[otu_to_lifehist,"range_gamete_km"])

par(mar = c(4,4,1,1))
stripchart(
  x = otu_sum_by_gamete_range, 
  method = "jitter", 
  pch = 21, 
  bg = rgb(1,1,1,alpha = 0.2), 
  las = 1
)

best_data_row <- max(which(!is.na(has_data)))
# THIS IS THE VALUE TO RETURN


life_history[,]
life_history[life_history$taxon_name == "Vetigastropoda",]
classification_by_query[as.character(classification_df$query[match(a_taxon, classification_df$taxname)])]


##########

for OTU (column)

from OTU TABLE

match on SEQ ID

from TAXON DF

print TAXON COLUMN

from TAXON DF

if TAXON

in TAXON COLUMN

from LIFE HISTORY TABLE

else

print LIFE HISTORY COLUMN

from LIFE HISTORY TABLE



match(
CLASSIFICATIONS_DF$TAXON
LIFE_HISTORY$TAXON
)

length(
match(
life_history$taxon_name,
classification_df$taxname
)
)

blast_df[taxa_vector, classification_col]
blast_df[, classification_col]

classification_df[, classification_col] %in%
life_history[, lh_taxon]


head(classification_df)
