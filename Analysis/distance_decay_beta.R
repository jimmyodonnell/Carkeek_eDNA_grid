# FIRST RUN 'load_data.r'

# dissimilarity by distance for reduced data set (mean abundance within sample for each OTU)
#-------------------------------------------------------------------------------

# REQUIRES:
# 1. dataframe called "metadata" with unique sequenced samples given in column "sample_id"
my_metadata <- metadata[["mean"]] # metadata[!duplicated(metadata[,colname_env_sample]),]
# 2. OTU table with rownames that correspond to aforementioned column "sample_id"
my_table <- otu_table[["mean"]] # mean, mean_unfilt, spvar, otu_named, as.binary(otu_mean), log, filt
rownames(my_table) # should be e.g. PCT-C-0000 etc, aka "env_sample_name"

EXPORT <- FALSE # export plots/files?

library(geosphere) # distm()
library(vegan) # vegdist()
library(propagate) # predictNLS
library(dplyr)
#-------------------------------------------------------------------------------
# calculate pairwise great circle distance between sampling locations using Haversine method
geo_dist <- as.dist(
  distm(
    x = my_metadata[,c(colname_lon, colname_lat)],
    fun = distHaversine
  )
)
attr(geo_dist, "Labels") <- my_metadata[, colname_env_sample]
# dimnames(geo_dist) <- list(my_metadata$env_sample_name, my_metadata$env_sample_name)

#-------------------------------------------------------------------------------
# calculate pairwise similarity of ecological communities
vegdist_method <- "bray"

distance_name <- switch(vegdist_method,
       bray     = "Bray_Curtis",
       morisita = "Morisita",
       horn     = "Morisita-Horn",
       jaccard  = "Jaccard",
       gower    = "Gower")

vegdist_methods <- c(
  "Bray_Curtis"   = "bray",
  "Morisita"      = "morisita",
  "Morisita_Horn" = "horn",
  "Jaccard"       = "jaccard",
  "Gower"         = "gower")

comm_sim <- list()
for(i in 1:length(vegdist_methods)){
  comm_sim[[names(vegdist_methods[i])]] <- vegdist(
    x = my_table, method = vegdist_methods[i], binary = FALSE)
}
for(i in 1:length(vegdist_methods)){
  comm_sim[[paste(names(vegdist_methods[i]), "bin", sep = "_")]] <- vegdist(
    x = my_table, method = vegdist_methods[i], binary = TRUE)
}

comm_sim <- lapply(comm_sim, function(x) 1 - x)

if(!(identical(attr(comm_sim[[1]], "Labels"), attr(geo_dist, "Labels")))){
	warning("Whoa there! the row/column names of the two distance matrices do not seem to add up. This is bad.")
}

# calculate replicate PCR similarity
otu_temp <- as.data.frame(prop(otu_table[["clean"]]))
PCR_similarities <- lapply(
  split(otu_temp, metadata[["clean"]][ , colname_env_sample]),
  function(x) {1 - vegdist(x, method = vegdist_method)}
)
rm(otu_temp)

#-------------------------------------------------------------------------------
# compare similarities of PCR replicates to environmental samples
PCR_mean <- round(mean(unlist(PCR_similarities)), digits = 3)
PCR_sd   <- round(  sd(unlist(PCR_similarities)), digits = 3)
env_mean <- round(mean(unlist(comm_sim["Bray_Curtis"])), digits = 3)
env_sd   <- round(  sd(unlist(comm_sim["Bray_Curtis"])), digits = 3)

text_similarity <- paste(
"PCR replicates within an environmental sample were extremely similar (",
PCR_mean, " plusminus ", PCR_sd,
") and far more similar than environmental samples (",
env_mean, " plusminus ", env_sd, ").",
sep = "")
par(mar = c(2,4,1,1))
boxplot(
list(PCR = unlist(PCR_similarities), environment = unlist(comm_sim["Bray_Curtis"])), 
las = 1, ylab = "Similarity"
)


#-------------------------------------------------------------------------------
# arrange the data for model fitting
model_data_full <- data.frame(
  dist2df(geo_dist),
  data.frame(lapply(comm_sim, as.vector))
)
# order the columns
model_data_full <- model_data_full[order(model_data_full[,"dist"]), ]


##########################
##########################
### --- STAN
##########################
##########################
library(rstan)

#-------------------------------------------------------------------------------
# MODEL IN STAN CODE
#-------------------------------------------------------------------------------
cat("
functions {#################################################################################
    real Beta(real mu , real phi, real y) {
      return( lgamma(phi) - lgamma(mu * phi) - lgamma(phi - phi*mu) + (mu*phi - 1) * log(y) + ((1-mu)*phi -1) * log(1-y))  ;
    }
}
data{
  // Define variables in data
    int<lower=0> N_obs ;
    vector<lower=0>[N_obs] bray;
    vector<lower=0>[N_obs] dist;
}

parameters {
// Define parameters to estimate
// intercept
real<lower=0,upper=1> beta_0;
// slope
real<lower=0> beta_1;
// precision
  real<lower=0> phi;
}

transformed parameters  {
  // Mean
  vector<lower=0,upper=1>[N_obs] mu;
  for(i in 1: N_obs){
    mu[i] = (beta_0 * beta_1) * pow(beta_1 + dist[i],-1) ;
  }
}

model {
// Priors

// Likelihood
  for(i in 1:N_obs){
    target += Beta(mu[i],phi,bray[i]) ;
  }
}
", file = "model.stan")

stan_data = list(
    "N_obs"     = nrow(model_data_full),
    "bray"      = model_data_full$Bray_Curtis,
    "dist"       = model_data_full$dist)

stan_pars = c(
  "beta_0",
  "beta_1",
  "mu",
  "phi")
  
Warm=1000
Iter=10000
N_CHAIN = 1
# Some options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
##############
stanMod = stan(file = 'model.stan',data = stan_data, 
               verbose = FALSE, chains = N_CHAIN, thin = 1, 
               warmup = Warm, iter = Warm + Iter, 
               #control = list(max_treedepth=10,adapt_delta=Adapt_delta,metric="diag_e"),
               pars = stan_pars)
               #sample_file = "/Users/ole.shelton/GitHub/Orca_Salmon/Output files/samp_file_mixed_test_E25_M2EST_noPUSO_vuln1_90FIN.csv",


pars <- extract(stanMod, permuted = TRUE)
traceplot(stanMod, pars = c("beta_0","beta_1","phi"))

DIST <- seq(0,5000,length.out=2000)

all <- NULL
for(i in 1:length(DIST)){  
  Y <- as.numeric(((pars$beta_0) * pars$beta_1) / ( pars$beta_1 + DIST[i]))
  temp <- c(Dist=DIST[i],Mean=mean(Y),SD=sd(Y),quantile(Y,probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975)))
  all <- rbind(all, temp)
}

all<- data.frame(all)

plot(Bray_Curtis~dist, data=model_data_full,ylim=c(0,1),xlim=c(0,4500))
par(new=T)
plot(Mean~Dist,data=all,ylim=c(0,1),xlim=c(0,4500),lwd=2,type="l",ylab="")
par(new=T)
plot(X97.5.~Dist,data=all,ylim=c(0,1),xlim=c(0,4500),lty=2,type="l",ylab="")
par(new=T)
plot(X2.5.~Dist,data=all,ylim=c(0,1),xlim=c(0,4500),lty=2,type="l",ylab="")












