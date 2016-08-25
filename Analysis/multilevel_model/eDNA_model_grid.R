#!/usr/bin/env Rscript

# This is code to estimate some model parameters for proportional abundance of DNA in a sample given reads in the file from a sequencer

library(R2jags)
library(reshape2) # acast: turn long-form dataframe into array

setwd("~/Documents/GoogleDrive/Kelly_Lab/Projects/Carkeek_eDNA_grid/Analysis")
setwd("/Users/jimmy.odonnell/Projects/Carkeek_eDNA_grid/Analysis")

data_dir <- file.path("..", "Data")
fig_dir <- file.path("..", "Figures")

# Let:
# i = taxon i (1:I); Carkeek grid dataset I = ? (limit to 10 for tests)
# j = PCR replicate j (1:J); Carkeek grid dataset J = 4
# k = position k from shore (1:K); Carkeek grid dataset K = 8
# l = position l along shore (1:L); L = 3

# Requires as input:
# 1. otu file
# 2. metadata file

# a table of counts of sequences (Z, length = ...)
# rows/rownames = samples, columns/colnames = taxa, cells = integer counts
counts_file_path <- file.path(data_dir, "OTUs_top10_4000m.csv")
counts_file_path <- file.path(data_dir, "otu_table_filtered_per_samp.csv")
counts_table <- read.csv(file = counts_file_path, row.names = 1)

# REMEMBER TO REMOVE CONTROL TAXON IF IT'S STILL THERE!
counts_table <- counts_table[,! names(counts_table) %in% "DUP_3"]


counts_table <- counts_table[,1:10] # use only the top ten for now

# if table is incorrectly oriented, transpose it:
# counts_table <- as.data.frame(t(counts_table))


# load metadata
metadata_file_path <- file.path(data_dir, "metadata_spatial.csv")
metadata <- read.csv(file = metadata_file_path, stringsAsFactors = FALSE)

# what is the name of the column containing the sample id (SEQUENCING samples)
colname_sampleid <- "sample_id"

# what is the name of the column containing PCR replicate levels
colname_pcr <- "PCR_replicate"

# what is the name of the column containing position in x direction (along shore) replicate levels
colname_posx <- "transect"

# what is the name of the column containing position in y direction (from shore) replicate levels
colname_posy <- "dist_from_shore"

# give a name to call the variable describing the things that were counted
varname_taxa <- "taxon"

# give a name to call the variable describing what the counts are
varname_counts <- "reads"


# of each taxon (length = I)
taxa <- colnames(counts_table)

# from each PCR replicate (length = J)
pcr <- as.character(unique(metadata[, colname_pcr]))

# at each position in x direction (length = K)
posx <- unique(metadata[, colname_posx])

# at each position in y direction (length = L)
posy <- as.character(unique(metadata[, colname_posy]))

# MAKE AN ARRAY OF DATA TO BE INPUT TO JAGS
# this really should be a function like this...
# JAGSarray3D <- function(table, split_vector){
#
# }

# *** NOTE *** #
# for the automated array approach to work, rows of metadata and otu table must line up!
if(nrow(counts_table) != nrow(metadata)){
	if( nrow(counts_table) > nrow(metadata) ){
		# extract and order the relevant rows of the OTU table
		counts_rel <- counts_table[metadata[, colname_sampleid], ]
		metadata_rel <- metadata
		# counts_table <- counts_rel
	} else if( nrow(counts_table) < nrow(metadata) ){
		# extract only the relevant rows of metadata
		metadata_rel <- metadata[match(rownames(counts_table), metadata[,colname_sampleid]),]
		counts_rel <- counts_table
		# metadata <- metadata_rel
	}
} else if(nrow(counts_table) == nrow(metadata)){
	counts_rel <- counts_table
	metadata_rel <- metadata
}

# match(metadata[, colname_sampleid], rownames(counts_table))

# make sure the order of the samples in the metadata and OTU table are the same
identical(
	rownames(counts_rel),
	metadata_rel[,colname_sampleid]
	)

# reshape into a long-format dataframe
data_l <- reshape(cbind(metadata_rel, counts_rel),
  varying = colnames(counts_rel), # aka taxa
  v.names = varname_counts,
  timevar = varname_taxa,
  times = colnames(counts_rel),
  new.row.names = 1:(nrow(counts_rel)*ncol(counts_rel)),
  direction = "long")[c(colname_sampleid, colname_pcr, colname_posx, colname_posy, varname_taxa, varname_counts)]

################################################################################
# THREE DIMENSIONAL ARRAY:

# this reorders the names of taxa
data_array <- acast(
				data = data_l,
				formula =  list(colname_pcr, varname_taxa, colname_posx),
				value.var = varname_counts
				)

# ... and thus this fucks up the names
# dimnames(data_array) <- list(
	# pcr = pcr,
	# taxa = taxa,
	# posx = posx
	# )

# so that if you're comparing this to the original way I created the array, this fails:
identical(counts, data_array)

# but this fixes it:
data_array <- data_array[,taxa,]
dimnames(data_array) <- list(
	pcr = pcr,
	taxa = taxa,
	posx = posx
	)
identical(counts, data_array)

data_array[ pcr[4] ,         ,         ]
data_array[        , taxa[8] ,         ]
data_array[        ,         , posx[2] ]

# # thus, you should be able to reference stuff like so:
counts[ pcr[4] ,         ,         ]
counts[        , taxa[8] ,         ]
counts[        ,         , posx[2] ]


################################################################################
# FOUR DIMENSIONAL ARRAY:

# this reorders the names of taxa
data_array <- acast(
				data = data_l,
				formula =  list(colname_pcr, varname_taxa, colname_posx, colname_posy),
				value.var = varname_counts
				)

# but this fixes it:
data_array <- data_array[,taxa,,]
dimnames(data_array) <- list(
	pcr = pcr,
	taxa = taxa,
	posx = posx,
	posy = posy
	)
identical(counts, data_array)

# # thus, you should be able to reference stuff like so:
data_array[ pcr[4] ,         ,         ,         ]
data_array[        , taxa[9] ,         ,         ]
data_array[        ,         , posx[3] ,         ]
data_array[        ,         ,         , posy[1] ]

counts <- data_array

# ORIGINAL APPROACH
# make a vector of the site for each sample name in the counts table
# posx_by_rowname <- metadata[ , colname_posx][
                        # match(rownames(counts_table),
                        # metadata[ , colname_sampleid])
                        # ]
# # unused for now
# # sampleid_by_posx <- split(x = metadata[, colname_sampleid], f = metadata[, colname_posx], drop = TRUE)
# counts_by_posx <- split(x = counts_table, f = posx_by_rowname)
# # create an array for JAGS called "counts"
# # where counts[j,i,k] is the counts in pcr replicate j of taxon i at site k
# counts <- array(
				# data = unlist(counts_by_posx),
				# dim = c(length(pcr), length(taxa), length(posx)),
				# dimnames = list(
				  # pcr = pcr,
				  # taxa = taxa,
				  # posx = posx
				# )
# )






# IMPORTANT
# JAGS won't do some basic R functions, so we need to give it the lengths explicitly
N_taxa <- length(taxa)
N_pcr <- length(pcr)
N_posx <- length(posx)
N_posy <- length(posy)


# lambda = mean of the Poisson dist that describes variation in counts
# beta = intercept (fixed, taxon-specific effect)

# eta = random effect of PCR for i,j,k,l
# sigma2 = variance of normal distribution describing eta (attributable to PCR)

# epsilon = random effect of position in X plane (along shore -- i.e. transect) for i,k,l
# tau2 = variance of normal distribution describing epsilon (attributable to position in X plane - k)

# delta = random effect of position in Y plane (distance from shore) i,l
# phi = variance of normal distribution describing delta (attributable to position in Y plane - l)



# pi = proportional abundance of each taxon

# model_loc <- "eDNA_model_grid.jags"


# Epsilon: i,k ; delta i,k ;
# jagsscript <- cat("
# jagsscript <- "
model {

    ## MODEL STRUCTURE
  for(l in 1:N_posy){
    for(k in 1:N_posx){
    	for(i in 1:N_taxa){
    	  for(j in 1:N_pcr){

            # Likelihood function
            # 3D: counts[j,i,k] ~ dpois(exp(lambda[j,i,k]))
    	      counts[j,i,k,l] ~ dpois(exp(lambda[j,i,k,l]))


            # GLM
            # 2D: lambda[i,j] <- beta_0 + beta[i] + eta[j,i]
            # alt format:
            # fixed[j,i] <- beta_0 + beta[i]
            # lambda[j,i] <- fixed[i] + eta[j,i]
    	    # 3D: lambda[j,i,k] <- beta_0 + beta[i] + eta[j,i,k] + epsilon[i,k]
    	    # 4D:
            lambda[j,i,k,l] <- beta_0 + beta[i] + eta[j,i,k,l] + epsilon[i,k,l] + delta[i,l]


            # random effect of PCR (for i,j)
            # note that precision = 1/variance and variance = sd^2
            # mu = mean, tau = precision
            # 2D: eta[j,i] ~ dnorm(0, 1/sigma2) # single site (k)
            # 3D: eta[j,i,k] ~ dnorm(0, 1/sigma2)
    	    # 4D:
            eta[j,i,k,l] ~ dnorm(0, 1/sigma2)


        }
      }
    }
  }
  for(l in 1:N_posy){
    for(k in 1:N_posx){
	    for(i in 1:N_taxa){
          # random effect of position in X plane (along shore -- i.e. transect) i,k,l
	      epsilon[i,k,l] ~ dnorm(0, 1/tau2[i])
	    }
    }
  }
  for(l in 1:N_posy){
    for(i in 1:N_taxa){
       # random effect of position in Y plane (distance from shore) i,l
      delta[i,l] ~ dnorm(0, 1/phi[i])
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
   	for(i in 1:N_taxa){
   		tau2[i] ~ dgamma(0.01, 0.01)
   	}

   	for(i in 1:N_taxa){
   	  phi[i] ~ dgamma(0.01, 0.01)
   	}
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


# ,
# file = model_loc)

jags_data <- list(
                "counts",
                "N_taxa",
                "N_pcr",
                "N_posx",
                "N_posy"
)

jags_params <- c(
                "beta",
                "sigma2",
                "tau2",
                "phi",
                "P",
                "beta_0"
)

# Set MCMC parameters
N_burn <- 0
N_iter <- 10000
N_chain <- 4

# run the MCMC
jags_out <- jags(
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
# attach.jags(jags_out, overwrite = TRUE)

# if you want to save:
# save(jags_out, file = "jagsoutput.RData")

# general plot
plot(jags_out)

# F this S. Keep getting the error "Error in coda.samples(jags_out, jags_params, N_iter = 1000) : attempt to apply non-function"
# coda.samples(jags_out, jags_params, N_iter = 1000)
# jags.samples(jags_out, jags_params, N_iter = 1000)











