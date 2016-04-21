# calculate dissimilarity
find_bad_PCR <- function(my_table, my_metadata, sample_id_column, grouping_column) {

	library(vegan) # vegdist
	# rownames must be unique sample IDs. get this from 1_data_prep.R

	# testing:
	# my_table <- otu_filt
	# my_metadata <- metadata
	# sample_id_column <- colname_sampleid
	# grouping_column <- colname_env_sample

	# At how many standard deviations above the mean dissimilarity
	# should a sample be excluded?
	THRESHOLD_SD <- 1
	if(!(identical(my_metadata[,sample_id_column], rownames(my_table)))){
		stop("hold on, the otu and metadata rows are not in the same order; disaster could ensue")
	}

	# vector of groups to which subsamples belong.
	# what is the name of the column containing the sample id from which subsamples were drawn (environmental samples)
	# colname_env_sample <- "env_sample_name"
	group_vector <- my_metadata[, grouping_column]


	dis_by_sample <- function(contingency_table, grouping_vector)
	{
	dis_by_group <- lapply(
	      # lapply(
	        split(as.data.frame(contingency_table), grouping_vector),
	        # matrix, ncol = ncol(contingency_table)
	      # ),
	      vegdist, method = "bray", upper = TRUE, diag = TRUE
	    )
	  return(dis_by_group)
	}

	my_dis <- dis_by_sample(my_table, group_vector)
	# unlist seems to conveniently remove 0s and redundant values for distance matrices
	my_dis_v <- unlist(my_dis)

	# calculate mean dissimilarity
	dis_mean <- mean(my_dis_v)
	dis_sd <- sd(my_dis_v)
	dissimilarity_cutoff <- dis_mean + (dis_sd * THRESHOLD_SD)

	# pdf(file = file.path(fig_dir, "BC_before_cleaning.pdf"))
	par(mar = c(4,7,1,1))
	stripchart(my_dis, method = "jitter", pch = 21, xlim = c(0,1),
	           xlab = "Pairwise Bray-Curtis Dissimilarity",
	           bg = rgb(0,0,0,alpha = 0.2), las = 1)


	# plot all
	# plot(my_dis_v, ylim = c(0,1),
		# ylab = "Within-sample dissimilarity (Bray-Curtis)",
		# xlab = "arbitrary sample index",
		# las = 1)
	line_colors <- c("black", "red")
	abline(v = c(dis_mean, dissimilarity_cutoff),
		lty = 2, col = line_colors)

	legend("topright", legend = c("mean", "threshold"),
		bty = "n", lty = 2, col = line_colors)

	# note high value for sample "lib_B_tag_GCGCTC" in PCT-C-0500
	# dev.off()

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

	# remove that sample?
	remove_bad_replicate <- TRUE
	if(remove_bad_replicate == TRUE){
	  metadata_cleaned <- my_metadata[my_metadata[,sample_id_column] != bad_replicate_name,]
	  # metadata_exp <- metadata_exp[metadata_exp$sample_id != bad_replicate_name,]
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
