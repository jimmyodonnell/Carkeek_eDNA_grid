# calculate dissimilarity
find_bad_replicate <- function(my_table, my_metadata,
	sample_id_column, grouping_column, threshold_sd = 1.5, 
    plot_results = TRUE) {

	if(!require(vegan)){stop("vegan package must be installed")}

	library(vegan) # vegdist
	# rownames must be unique sample IDs. get this from 1_data_prep.R

	# convert the data matrix to proportional abundance within sample
	prop_table <- my_table/rowSums(my_table)

	# testing:
	# my_table <- otu_filt
	# my_metadata <- metadata
	# sample_id_column <- colname_sampleid
	# grouping_column <- colname_env_sample

	# threshold_sd: At how many standard deviations above the mean dissimilarity
	# should a sample be excluded?

	if(!(identical(my_metadata[,sample_id_column], rownames(my_table)))){
		stop("hold on, the otu and metadata rows are not in the same order; disaster could ensue")
	}

	# vector of groups to which subsamples belong.
	# what is the name of the column containing the sample id from which subsamples were drawn (environmental samples)
	# colname_env_sample <- "env_sample_name"
	group_vector <- my_metadata[, grouping_column]


	my_dis <- lapply(
	  split(as.data.frame(prop_table), group_vector),
	  vegdist, method = "bray", binary = FALSE, diag = TRUE, upper = TRUE
	)

	# unlist seems to conveniently remove 0s and redundant values for distance matrices
	my_dis_v <- unlist(my_dis)

	# calculate mean dissimilarity
	dis_mean       <- mean(my_dis_v)
	dis_mean_round <- round(dis_mean, digits = 3)
	dis_sd         <- sd(my_dis_v)
	dis_sd_round   <- round(dis_sd, digits = 3)
	print(
		paste(
			"dissimilarity among replicates within sample (mean, sd):", 
			paste(
				dis_mean_round, 
				dis_sd_round, 
				collapse = " ")))
	dissimilarity_cutoff <- dis_mean + (dis_sd * threshold_sd)
    dissimilarity_cutoff_round <- round(dissimilarity_cutoff, digits = 3)

	# which sample group contains dissimilarities over the threshold?
	bad_env_sample <- which(
		sapply(
			lapply(my_dis,
			function(x) which(x > dissimilarity_cutoff)
			), length) > 0)


	# Which samples have high dissimilarity among other samples from the same environmental sample?
	# what are the names of those bad comparisons within that sample
	dis_names <- rownames(which(as.matrix(my_dis[[bad_env_sample]]) > dissimilarity_cutoff, arr.ind = TRUE))

	# what is the name of the replicate occuring more than once?
	bad_replicate_name <- names(table(dis_names))[table(dis_names) > 1]

	if(plot_results){
		color_vector <- rep(1, length(my_dis))
		color_vector[bad_env_sample] <- "white"
		par(mar = c(4,7,1,1))
		stripchart(my_dis, method = "jitter", pch = 21, xlim = c(0,1),
		           xlab = "Pairwise Bray-Curtis Dissimilarity",
		           bg = rgb(0,0,0,alpha = 0.2), #col = color_vector,
		           las = 1)
	    abline(h = seq(from = 0.5, to = 60.5, by = 1), col = "gray")
	
		# attempting to color the points red over the threshold
		# bad_x <- as.vector(my_dis[[bad_env_sample]])
		# bad_y <- rep(bad_env_sample, length(temp))
		# points(jitter(bad_x), jitter(bad_y), col = (as.numeric(bad_x > dissimilarity_cutoff) + 1))
	
		# plot all
		# plot(my_dis_v, ylim = c(0,1),
			# ylab = "Within-sample dissimilarity (Bray-Curtis)",
			# xlab = "arbitrary sample index",
			# las = 1)
		line_colors <- c("black", "red")
		abline(v = c(dis_mean, dissimilarity_cutoff),
			lty = 2, col = line_colors)
	
		legend("topright", 
			legend = c(
				paste("mean =", dis_mean_round), 
				paste("threshold =", dissimilarity_cutoff_round)),
			bty = "n", lty = 2, col = line_colors)
	}

	# note high value for sample "lib_B_tag_GCGCTC" in PCT-C-0500

	# remove the bad sample?
	remove_bad_replicate <- TRUE
	if(remove_bad_replicate == TRUE){
	  metadata_cleaned <- my_metadata[my_metadata[,sample_id_column] != bad_replicate_name,]
	  my_table_cleaned <- my_table[rownames(my_table) != bad_replicate_name,]
	  print(paste("removed outlier PCR:", bad_replicate_name))
	}

	return(list(my_table_cleaned, metadata_cleaned))

}

# my_table <- my_table_cleaned
# if(!(identical(metadata[,colname_sampleid], rownames(my_table)))){
	# warning("hold on, the otu and metadata rows are not in the same order; disaster could ensue")
# }
# vector of groups to which samples belong.
# group_vector <- metadata$env_sample_name
# my_dis <- dis_by_sample(my_table, group_vector)

# pdf(file = file.path(fig_dir, "BC_after_cleaning.pdf"))
# par(mar = c(4,7,1,1))
# stripchart(my_dis, method = "jitter", pch = 21, las = 1,
           # xlab = "Pairwise Bray-Curtis Dissimilarity", bg = rgb(0,0,0,alpha = 0.2))
# dev.off()

# calculate mean dissimilarity
# dis_mean <- lapply(my_dis, function(x) rowMeans(as.matrix(x)))
# max(sapply(dis_mean, mean))

####################################################################################
# ALTERNATE
# Use on mean OTU abundance per water sample (i.e. don't compare among PCR replicates)
# vegdist(otu_mean, method = "bray", binary = FALSE)
