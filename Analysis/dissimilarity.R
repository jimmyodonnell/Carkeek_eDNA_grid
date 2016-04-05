# calculate dissimilarity

library(vegan) # vegdist
# rownames must be unique sample IDs. get this from 1_data_prep.R
otu_table <- otu_filt


# vector of groups to which samples belong.
group_vector <- metadata$env_sample_name

dis_by_sample <- function(contingency_table, grouping_vector)
{
  return(
    lapply(
      lapply(
        split(otu_table, group_vector), 
        matrix, ncol = ncol(otu_table)
      ), 
      vegdist, method = "bray", upper = TRUE, diag = TRUE
    )
  )
}

my_dis <- dis_by_sample(otu_table, group_vector)



par(mar = c(4,6,1,1))
stripchart(my_dis, method = "jitter", pch = 21, las = 1, 
           xlab = "Pairwise Bray-Curtis Dissimilarity", bg = rgb(0,0,0,alpha = 0.2))

# note high value for sample "lib_B_tag_GCGCTC" in PCT-C-0500

# my_dis[1]
# lapply(my_dis, function(x) max())