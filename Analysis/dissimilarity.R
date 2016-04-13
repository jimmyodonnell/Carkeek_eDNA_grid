# calculate dissimilarity

library(vegan) # vegdist
# rownames must be unique sample IDs. get this from 1_data_prep.R
my_table <- otu_filt

if(!(identical(metadata[,colname_sampleid], rownames(my_table)))){
	warning("hold on, the otu and metadata rows are not in the same order; disaster could ensue")
}

# vector of groups to which subsamples belong.
# what is the name of the column containing the sample id from which subsamples were drawn (environmental samples)
colname_env_sample <- "env_sample_name"
group_vector <- metadata[, colname_env_sample]


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

# which sample has the highest BC dissimilarity?
dis_high <- which.max(sapply(my_dis, max))
as.dist(my_dis[[dis_high]], upper = FALSE, diag = FALSE)
mean(as.dist(my_dis[[dis_high]], upper = FALSE, diag = FALSE)[1:3])

# calculate mean dissimilarity
dis_mean <- lapply(my_dis, function(x) rowMeans(as.matrix(x)))

# Which samples have high dissimilarity among other samples from the same environmental sample?
dis_mean_v <- unlist(dis_mean)
names(dis_mean_v) <- do.call(c, lapply(dis_mean, names))

# One "bad" replicate will have a high mean BC dissimilarity to all other replicates, but it will also elevate the BC dissimilarity of the other replicates
round(sort(dis_mean_v), digits = 3)
plot(dis_mean_v)
# pdf(file = file.path(fig_dir, "BC_before_cleaning.pdf"))
par(mar = c(4,7,1,1))
stripchart(my_dis, method = "jitter", pch = 21, las = 1, 
           xlab = "Pairwise Bray-Curtis Dissimilarity", bg = rgb(0,0,0,alpha = 0.2))
# dev.off()

# note high value for sample "lib_B_tag_GCGCTC" in PCT-C-0500
# in this case, 0.2 seems to be a reasonable cutoff beyond which we think a PCR replicate is "bad"
names(which(dis_mean_v > 0.2))
names(which.max(dis_mean_v))

# remove that sample?
remove_outlier <- TRUE
if(remove_outlier == TRUE){
  outlier <- names(which(dis_mean_v > 0.2))
  metadata <- metadata[metadata$sample_id != outlier,]
  metadata_exp <- metadata_exp[metadata_exp$sample_id != outlier,]
  otu_filt <- otu_filt[rownames(otu_filt) != outlier,]
  print("removed outlier PCR")
}

my_table <- otu_filt
if(!(identical(metadata[,colname_sampleid], rownames(my_table)))){
	warning("hold on, the otu and metadata rows are not in the same order; disaster could ensue")
}
# vector of groups to which samples belong.
group_vector <- metadata$env_sample_name
my_dis <- dis_by_sample(my_table, group_vector)


# my_dis[1]
# lapply(my_dis, function(x) max())

# pdf(file = file.path(fig_dir, "BC_after_cleaning.pdf"))
par(mar = c(4,7,1,1))
stripchart(my_dis, method = "jitter", pch = 21, las = 1, 
           xlab = "Pairwise Bray-Curtis Dissimilarity", bg = rgb(0,0,0,alpha = 0.2))
# dev.off()

# calculate mean dissimilarity
dis_mean <- lapply(my_dis, function(x) rowMeans(as.matrix(x)))
max(sapply(dis_mean, mean))

####################################################################################
# ALTERNATE
# Use on mean OTU abundance per water sample (i.e. don't compare among PCR replicates)
vegdist(otu_mean, method = "bray", binary = FALSE)