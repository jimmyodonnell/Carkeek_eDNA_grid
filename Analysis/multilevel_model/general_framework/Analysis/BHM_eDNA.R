# Estimation of some underlying/unobserved quantity from some other observed quantity
# (e.g. relative abundance from eDNA counts)

# Set this to the analysis directory of this project
analysis_dir <- "~/Documents/GoogleDrive/Kelly_Lab/Projects/Carkeek_eDNA_grid/Analysis"



setwd(analysis_dir)
data_dir <- file.path("..", "Data")
fig_dir <- file.path("..", "Figures")


############################################################################
# The following requires as input:
# 1. A matrix of counts (integers) of DNA sequences, where:
# 		rows = samples
# 		columns = taxa whose names are given in column names
# eg: dput(mydata) 
# structure(c(18094L, 18136L, 21198L, 13643L, 14198L, 12219L, 12703L, 
# 13873L, 11773L, 12022L, 11865L, 13554L, 13840L, 9871L, 10999L, 
# 10070L, 10827L, 8920L, 9708L, 10033L, 9303L, 6802L, 7853L, 7927L, 
# 3118L, 3344L, 4107L), .Dim = c(3L, 9L), .Dimnames = list(NULL, 
    # c("taxon_A", "taxon_B", "taxon_C", "taxon_D", "taxon_E", 
    # "taxon_F", "taxon_G", "taxon_H", "taxon_I")))


# If loading from a file
# counts <- as.matrix(read.csv(file = file.path(data_dir, "counts.csv")))
# taxon_names <- colnames(counts)

# If reading from a currently existing object (e.g. you simulated some eDNA data)
# counts <- as.matrix(mydata)

# make into a vector
counts_vec <- as.vector(counts)
names(counts_vec) <- rep(x = colnames(counts), each = nrow(counts))




############
# Create Design Matrix for the categorical effects
############

# Should have nrow = length of the input 

cat_mat <- matrix(data = 0, nrow = length(counts_vec), ncol = ncol(counts))

# Not sure why, but you set the first species to equal 1...
# I think this might be because everything else is expressed relative to this?
cat_mat[,1]	<- 1

# Changed Ole's original version because it will break if you don't use the "right" taxon names (i.e. with an integer corresponding to their index in the third position == A.3)
for(i in 1:length(taxon_names)){
	same_taxon <- which(rep(1:length(taxon_names), each = N_samples_small) == i)
	cat_mat[same_taxon,i]	<-	1
	# Ole's solution is more readable (x <- 1), but you can also make use of the fact that logical vectors (T/F) can be expressed as 1/0 using as.integer
	# as.integer(rep(1:length(taxon_names), each = N_samples_small) == i)
}
colnames(cat_mat)	<-	taxon_names




############
# Create Design Matrix for the variance
############

# Different Variance for each species
var_mat	<-	matrix(data = 0, nrow = length(counts_vec), ncol = length(taxon_names))

for(i in 1:length(taxon_names)){
	same_taxon <- which(rep(1:length(taxon_names), each = nrow(counts)) == i)
	var_mat[same_taxon, i]	<-	1
}
colnames(var_mat) <- taxon_names

######
# Establish beta parameter vector
######

betas	<- rep(0, ncol(cat_mat))


##########################################################################################
##### Run model WITHOUT random effect
##########################################################################################


N_cov	<-	length(taxon_names)
N_obs	<-	length(counts_vec)

# Identity matrix for all the observations
iden_mat <- diag(x = 1, nrow = N_obs)


# Let's see if I can put the jags variables down before writing the JAGS script

# WHAT IS y? (counts_vec)
jags_data <- list("counts_vec", "N_obs", "N_cov", "cat_mat")

# WHAT IS P?
jags_params <- c("betas", "P")

# I think this is telling where the script lives?
model_loc <- c("my_jags_script.txt")

# Set the number of iterations to discard before storing them
N_burn <- 10000

# Set the total number of iterations to conduct (total? or after burnin?)
N_iter <- 10000


# write the JAGS script:
jagsscript <- cat(
"model {
	for(i in 1:N_obs){
		counts_vec[i] ~ dpois(exp(lambda[i]))
	}
	
	lambda <- cat_mat %*% betas #+ eta
	
	### Derived Quantities
	p[1] <- exp(betas[1])
	for(i in 2:N_cov){
		p[i] <- exp(betas[1] + betas[i])
	}
	
	# This can probably just be: P <- p / sum(p)
	for(i in 1:N_cov){
		P[i] <- p[i] / sum(p)
	}
	
	### PRIORS
	for(j in 1:N_cov){
		betas[j] ~ dnorm(0, 0.001)
	}
}", file = model_loc)

#) # cat doesn't close until I make this parentheses...
#)

# What is this doing?
my_jags_model <- jags(
					data = jags_data, 
					inits = NULL, 
					parameters.to.save = jags_params, 
					model.file = model_loc, 
					n.chains = 3, 
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


# What does this do?
attach.jags(my_jags_model, overwrite = TRUE)

























##########################################################################################
##### Run model WITH random effect
##########################################################################################

