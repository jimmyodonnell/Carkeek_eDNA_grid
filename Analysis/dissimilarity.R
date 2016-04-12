# calculate dissimilarity

library(vegan) # vegdist
# rownames must be unique sample IDs. get this from 1_data_prep.R
my_table <- otu_filt

if(!(identical(metadata[,colname_sampleid], rownames(my_table)))){
	warning("hold on, the otu and metadata rows are not in the same order; disaster could ensue")
}
# vector of groups to which samples belong.
group_vector <- metadata$env_sample_name


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

# calculate mean dissimilarity
dis_mean <- lapply(my_dis, function(x) rowMeans(as.matrix(x)))

# Which samples have high dissimilarity among other samples from the same environmental sample?
dis_mean_v <- unlist(dis_mean)
names(dis_mean_v) <- do.call(c, lapply(dis_mean, names))
names(which.max(dis_mean_v))

par(mar = c(4,6,1,1))
stripchart(my_dis, method = "jitter", pch = 21, las = 1, 
           xlab = "Pairwise Bray-Curtis Dissimilarity", bg = rgb(0,0,0,alpha = 0.2))

# note high value for sample "lib_B_tag_GCGCTC" in PCT-C-0500

# my_dis[1]
# lapply(my_dis, function(x) max())


####################################################################################
# ALTERNATE
# Use on mean OTU abundance per water sample (i.e. don't compare among PCR replicates)
vegdist(otu_mean, method = "bray", binary = FALSE)