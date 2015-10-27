#!/usr/bin/env Rscript

# This is code to estimate some model parameters for proportional abundance of DNA in a sample given reads in the file from a sequencer

library(R2jags)

setwd("~/Documents/GoogleDrive/Kelly_Lab/Projects/Carkeek_eDNA_grid/Analysis")

data_dir <- file.path("..", "Data")
fig_dir <- file.path("..", "Figures")

# Let:
# i = taxon i (1:I); Carkeek grid dataset I = ? (limit to 10 for tests)
# j = PCR replicate j (1:J); Carkeek grid dataset J = 4
# k = location k (1:K); Carkeek grid dataset K = 24


# Requires as input:
# 1. otu file
# 2. metadata file

# a table of counts of sequences (Z, length = ...)
counts_file_path <- file.path(data_dir, "OTUs_top10_4000m.csv")
counts_table <- read.csv(file = counts_file_path, row.names = 1)
# rows/rownames = samples, columns/colnames = taxa, cells = integer counts
# if table is incorrectly oriented, transpose it:
# counts_table <- as.data.frame(t(counts_table))

# load metadata
metadata_file_path <- file.path(data_dir, "metadata_spatial.csv")
metadata <- read.csv(file = metadata_file_path, stringsAsFactors = FALSE)

# what is the name of the column containing the sample id (SEQUENCING samples)
colname_sampleid <- "sample_id"

# what is the name of the column containing PCR replicate levels
colname_pcr <- "PCR_replicate"

# what is the name of the column containing relevant site replicate levels
colname_site <- "transect"

# of each taxon (length = I)
taxa <- colnames(counts_table)

# from each PCR replicate (length = J)
pcr <- unique(metadata[, colname_pcr])

# at each site (length = K)
sites <- unique(metadata[, colname_site])








# JAGS won't do some basic R functions, so we need to give it the lengths explicitly
N_taxa <- length(taxa)
N_pcr <- length(pcr)
N_sites <- length(sites)


# make a vector of the site for each sample name in the counts table
site_by_rowname <- metadata[ , colname_site][
  match(rownames(counts_table), metadata[ , colname_sampleid])
]

# unused for now
# sampleid_by_site <- split(x = metadata[, colname_sampleid], f = metadata[, colname_site], drop = TRUE)

counts_by_site <- split(x = counts_table, f = site_by_rowname)






# create an array for JAGS called "counts"
# where counts[j,i,k] is the counts in pcr replicate j of taxon i at site k
counts <- array(
				data = unlist(counts_by_site),
				dim = c(length(pcr), length(taxa), length(sites)),
				dimnames = list(
				  pcr = pcr,
				  taxa = taxa,
				  sites = sites
				)
)

# this really should be a function like this...
# JAGSarray3D <- function(table, split_vector){
#
# }

# thus, you should be able to reference stuff like so:
counts[ pcr[4] ,         ,          ]
counts[        , taxa[8] ,          ]
counts[        ,         , sites[2] ]

# lambda = mean of the Poisson dist that describes variation in counts
# beta = intercept
# eta = random effect for i,j,k
# sigma2 = variance of normal distribution describing eta (attributable to PCR)

# epsilon = random effect for i,k
# tau2 = variance of normal distribution describing epsilon (attributable to location - k)

# pi = proportional abundance of each taxon

model_loc <- "eDNA_model_grid.jags"

# jagsscript <- cat("
jagsscript <- "
model {

    ## MODEL STRUCTURE
  for(k in 1:N_sites){
    	for(i in 1:N_taxa){
    	  for(j in 1:N_pcr){

            # Likelihood function
            counts[j,i,k] ~ dpois(exp(lambda[j,i,k]))

            # GLM
            lambda[j,i,k] <- beta_0 + beta[i] + eta[j,i,k] + epsilon[j,i,k]
            # single site: lambda[i,j] <- beta_0 + beta[i] + eta[j,i]
            # alt format:
            # fixed[j,i] <- beta_0 + beta[i]
            # lambda[j,i] <- fixed[i] + eta[j,i]

            # random effect for i,j
            # note that precision = 1/variance and variance = sd^2
            # mu = mean, tau = precision
            eta[j,i,k] ~ dnorm(0, 1/sigma2)
            # eta[j,i] ~ dnorm(0, 1/sigma2) # single site (k)

            # random effect for i,k
            epsilon[j,i,k] ~ dnorm(0, 1/tau2)

        }
      }
    }

    ## PRIORS
    # *NOTE: in JAGS gamma, shape = r and rate = mu (lambda?)

    # the general intercept (mean of counts of all taxa for all PCR replicates)
    beta_0 ~ dnorm(0, 1/1000)

    # beta[1] must be zero, because for taxa[1]
    beta[1] <- 0

    for(i in 2:N_taxa){
        # mu = mean, tau = precision
        beta[i] ~ dnorm(0, 1/1000)
    }

    # to allow variance attributable to PCR to vary among species,
    # can make this species-specific by adding [i] and looping
    # for(i in 1:N_taxa){
    # in JAGS gamma, shape = r and rate = mu (lambda?)
   	sigma2 ~ dgamma(0.01, 0.01)

   	# in JAGS gamma, shape = r and rate = mu (lambda?)
   	tau2 ~ dgamma(0.01, 0.01)

   ## DERIVED QUANTITIES

    # multinomial poisson transformation
    # i.e. estimated proportion of taxa[i]
    for(i in 1:N_taxa){
        p[i] <- exp(beta_0 + beta[i])
    }

    for(i in 1:N_taxa){
        P[i] <- p[i] / sum(p)
    }

}
"
"


# ,
# file = model_loc)

jags_data <- list(
                "counts",
                "N_taxa",
                "N_pcr",
                "N_sites"
)

jags_params <- c(
                "beta",
                "sigma2",
                "tau2",
                "P",
                "beta_0"
)

# Set MCMC parameters
N_burn <- 0
N_iter <- 100000
N_chain <- 3

# run the MCMC
my_jags <- jags(
				data = jags_data,
				inits = NULL,
				parameters.to.save = jags_params,
				model.file = textConnection(jagsscript),
				# model.file = model_loc,
				n.chains = N_chain,
				n.iter = N_iter + N_burn,
				n.burnin = N_burn, #floor(n.iter/2)
				n.thin = 1, # max(1, floor((n.iter - n.burnin)/1000))
				DIC = TRUE,
				working.directory = NULL,
				jags.seed = 123,
				refresh = N_iter/50,
				progress.bar = "text",
				digits = 5,
				RNGname = c("Wichmann-Hill", "Marsaglia-Multicarry", "Super-Duper", "Mersenne-Twister"),
				jags.module = c("glm","dic")
)

# Attach the jags output
attach.jags(my_jags, overwrite = TRUE)

# if you want to save:
# save(my_jags, file = "jagsoutput.RData")

# general plot
plot(my_jags)

# F this S. Keep getting the error "Error in coda.samples(my_jags, jags_params, N_iter = 1000) : attempt to apply non-function"
# coda.samples(my_jags, jags_params, N_iter = 1000)
# jags.samples(my_jags, jags_params, N_iter = 1000)

str(my_jags)
names(my_jags) # "model" "BUGSoutput" "parameters.to.save" "model.file" "n.iter" "DIC"

str(my_jags$BUGSoutput)
names(my_jags$BUGSoutput)

# simulations list
names(my_jags$BUGSoutput$sims.list)

# simulations array
str(my_jags$BUGSoutput$sims.array)


dim(my_jags$BUGSoutput$sims.array)
dimnames(my_jags$BUGSoutput$sims.array)
dimnames(my_jags$BUGSoutput$sims.array)[[3]] # parameter names

dim(my_jags$BUGSoutput$sims.array[,,])




mcmc_summary <- my_jags$BUGSoutput$summary








################################################################################################
# PLOT TRACES
################################################################################################

mcmc_array <- my_jags$BUGSoutput$sims.array

dim(mcmc_array)

mcmc_iter <- dim(mcmc_array)[1]
mcmc_chains <- dim(mcmc_array)[3]
mcmc_params <- dimnames(mcmc_array)[[3]]
thin_factor <- 10

trace_layout <- layout(

		mat = matrix(
					data = 1:length(mcmc_params), 
					nrow = 6, 
					byrow = TRUE
					)
	)

layout.show(trace_layout)

for(param in mcmc_params){

	mcmc_samples <- mcmc_array[,,param]

	plot(
		x= 0, 
		pch = '', 
		xlim = c(0, nrow(mcmc_samples)/thin_factor), 
		ylim = c(min(mcmc_samples), max(mcmc_samples)), 
		xlab = "MCMC iteration", 
		ylab = "value", 
		main = param
		)
		
	for(i in 1:ncol(mcmc_samples)){
	
		thinned_chain <- seq(from = 1, to = nrow(mcmc_samples), by = thin_factor)

		points(
			x = mcmc_samples[thinned_chain,i], 
			type = "l", 
			col = rgb(0,0,0, alpha = 0.3)
			)

	}

}


# OTHER PLOTS
names(my_jags$BUGSoutput$sims.list)

dim(my_jags$BUGSoutput$sims.list$deviance) # N_iter*N_chain by 1


boxplot(
	x = my_jags$BUGSoutput$sims.list$deviance,
	main = "Deviance"
	)

dim(my_jags$BUGSoutput$sims.list[[1]]) # matrix of dim 300000 (N_iter*N_chain) by 10 (N_taxa)

# so, to plot the estimated proportion of taxon 1 from all samples:
boxplot(my_jags$BUGSoutput$sims.list[[1]][,1])

# or, plot the estimates of P for all species:

pdf(file = file.path(fig_dir, "P_boxplot.pdf"))
boxplot(
	x = my_jags$BUGSoutput$sims.list[[1]], 
	main = "Posterior estimates of P", 
	xlab = "Taxa", 
	xaxt = "n", 
	# names = taxa, 
	las = 1, 
	ylim = c(0,1), 
	# cex.axis = 0.5
	)
axis(side = 1, at = 1:length(taxa), labels = taxa, cex.axis = 0.7)
dev.off()
dim()
