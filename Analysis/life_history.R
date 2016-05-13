# add life history data

life_history_file <- "life_history.csv"
life_history <- read.csv(
  file = file.path(data_dir, life_history_file),
  stringsAsFactors = FALSE
)

head(life_history)
lh_taxon <- "taxon_name"

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
# their index in the life history data
classification_df$lifehist_row <- match(
	classification_df$taxname, life_history$taxon_name
)
head(classification_df)


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
