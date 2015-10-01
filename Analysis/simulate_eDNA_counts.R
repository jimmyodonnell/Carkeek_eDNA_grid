############################################################################
# This script simulates some data that resembles eDNA sequence data
############################################################################

# Original code by Ole Shelton
# Hacked up by Jimmy O'Donnell
# Comments mostly reflect Jimmy's imperfect understanding of the underlying concepts.

library(gtools) # gtools::rdirichlet
library(R2jags) # R2jags::jags



# Suppose we have an aquarium which contains N_tax taxa
N_tax <- 9



# Give them names
taxon_names <- paste("taxon_", LETTERS[1:N_tax], sep = "")


# The distribution of DNA molecules derived from each of the taxa might look like the following
# Ole originally called these alpha because they eventually become the alpha parameter of the Dirichlet distribution, which sorta determines how likely each component is to be drawn
# I call these "true" because they represent the "true" abundance of DNA from each taxon in the sample.
# Set seed to be able to reproduce draws from a probability distribution (gets reset each time one of these functions is used)
# Instead of runif, you might try using rexp or random from the Pareto distribution to more accurately resemble "real" eDNA results
set.seed(407744)
DNA_truth		<- sort(round(runif(N_tax,1,100)), decreasing = TRUE)
names(DNA_truth) <- taxon_names

# those same values expressed as proportions
DNA_truth_prop <- DNA_truth/sum(DNA_truth)





############################################################################
# The actual simulation 
############################################################################

# This simulates data that might be expected from high throughput sequencing of PCR amplicons generated from environmental samples: proportional abundance of sequences from each of the taxa/OTUs/dups
# draw samples from the Dirichlet distribution, using alpha given by DNA_truth to approximate the sampling of DNA molecules from the environment
# Set seed to be able to reproduce draws from a probability distribution (gets reset each time one of these functions is used)


# If you were to take N_samples samples, what kind of answer would you expect?
N_samples_small	<-	3
set.seed(407744)
eDNA_counts_small <- rdirichlet(N_samples_small, DNA_truth)
colnames(eDNA_counts_small) <- taxon_names


############
# OR...
############
N_samples_big		<- 100000
set.seed(407744)
eDNA_counts_big <- rdirichlet(N_samples_big, DNA_truth)
colnames(eDNA_counts_big) <- taxon_names


# Set whichever of these as the plot input:
plot_input <- eDNA_counts_big

# Create plot layout (this looks complicated to allow it to be flexible to different taxa numbers)
plot_cols <- ceiling(sqrt(ncol(plot_input)))
plot_rows <- ceiling(ncol(plot_input)/plot_cols)
plot_layout <- matrix(data = 1:(plot_cols*plot_rows), nrow = plot_rows, byrow = TRUE)
layout(plot_layout)

# plot histograms of the frequency of proportions sampled for each taxon.
for (i in 1:ncol(plot_input)){
  hist(plot_input[, i], main = taxon_names[i], xlab = "proportion")
  abline(v = DNA_truth_prop[i], col = "red")
}


# In a similar way, we can use each of the samples to inform the probabilities of draws from a multinomial distribution. This produces data that is functionally indistinguishable from that presented above.
# Using each row of the "true proportion" data frame as probabilities...
# essentially: take a draw of 'N_draws' marbles, 
# and put them into each of some number of bins 
# with probability of going into each bin given by the "true proportion" row
# (repeat n times)

N_draws		<-	100000

set.seed(407744)
counts_mat <- apply(
					X = eDNA_counts_small, 
					MARGIN = 1, 
					FUN = function(x) 
						rmultinom(size = N_draws, prob = x, n = 1)
				)

# assign names
rownames(counts_mat) <- taxon_names

# transpose
counts <- t(counts_mat)

# want to save the file?
write.csv(x = counts, file = "counts.csv", row.names = FALSE, quote = FALSE)

# or to prepare for next step (MCMC),
mydata <- counts

############################################################################
# This marks the end of generating what I think of as typical eDNA sequence data
############################################################################















































