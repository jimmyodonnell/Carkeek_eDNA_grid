# calculate dissimilarity

library(vegan) # vegdist
# rownames must be unique sample IDs. get this from 1_data_prep.R
my_table <- otu_filt

vegdist(otu_filt)
# vector of groups to which samples belong.
group_vector <- metadata$env_sample_name

dis_by_sample <- function(contingency_table, grouping_vector)
{
  return(
    lapply(
      lapply(
        split(contingency_table, grouping_vector), 
        matrix, ncol = ncol(contingency_table)
      ), 
      vegdist, method = "bray", upper = TRUE, diag = TRUE
    )
  )
}

my_dis <- dis_by_sample(my_table, group_vector)

which.max(rowSums(do.call(rbind, lapply(my_dis, function(x) rowMeans(as.matrix(x))))))
which.max(do.call(rbind, lapply(my_dis, function(x) rowMeans(as.matrix(x)))))

metadata[,"sample_id"][5]
my_dis[[1]]

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