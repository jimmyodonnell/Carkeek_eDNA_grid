
classifications_file <- "classifications20160329.RData"

# I don't know how this works, I am just clicking to load the data
# classifications <- load(file.path(data_dir, classifications_file))

classifications


the_ranks <- lapply(lapply(classifications, function(x) x[2]), "[[", 1)
the_names <- lapply(classifications, "[[", 1)

# unique(unlist(lapply(lapply(classifications, function(x) x[2]), "[[", 1), use.names = FALSE))

classification_df <- data.frame(
         query = rep(names(the_names), lapply(the_names, length)),
         rankname = unlist(the_ranks), 
         taxname = unlist(the_names), 
         row.names = NULL
         )

current_time <- gsub(pattern = " ", replacement = "_", x = Sys.time())
write.csv(classification_df,file = file.path(data_dir, paste("classification_", current_time, ".csv", sep = "")))

#----------------------------------------------------------------------------------------
# Get names of taxonomic ranks (e.g. "kingdom", "subphylum", etc)
#----------------------------------------------------------------------------------------
ranknames <- getranknames()
unique_ranks <- sort(as.numeric(unique(ranknames[,"rankid"])))
all_ranks <- tolower(ranknames[match(as.character(unique_ranks), ranknames[,"rankid"]),"rankname"]) # c(, "no rank")
all_ranks_full <- c(rbind(all_ranks, paste("below-", all_ranks, sep = "")))



