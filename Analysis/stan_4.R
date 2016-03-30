data_dir <- file.path("..", "Data")

data_file <- file.path(data_dir, "data_full_long_20.csv")

data_full <- read.csv(data_file)

mydata <- data_full

colname_level1 <- "sequenced_sample"
colname_level2 <- "transect_line"
colname_level3 <- "transect_position"
colname_level4 <- "OTU"
colname_predictor <- "dist_from_shore"
colname_outcome <- "count"

lookup_3in2 <- unique(mydata[,c(colname_level2, colname_level3)])[,colname_level3]
lookup_4in3 <- as.character(unique(mydata[,c(colname_level3, colname_level4)])[,colname_level4])

# plot(mydata[,colname_predictor], mydata[,colname_outcome], xlab = "distance from shore", ylab = "number of reads", las = 1)

stan_data <- list(Ni          = length(unique(mydata[,colname_level1])), 
                  Nj          = length(unique(mydata[,colname_level2])), 
                  Nk          = length(unique(mydata[,colname_level3])), 
                  Nl          = length(unique(mydata[,colname_level4])), 
                  level2      = mydata[,colname_level2],
                  level3      = mydata[,colname_level3],
                  level4      = mydata[,colname_level4],
                  lev3ForLev2 = lookup_3in2,
                  lev4ForLev2 = lookup_4in3,
                  Y_ijk       = mydata[,colname_outcome],
                  X_1ijk      = mydata[,colname_predictor]
                  )

lm_out <- lm(mydata[,colname_outcome] ~ mydata[,colname_predictor] )
summary(lm_out)

#-------------------------------------------------------------------------------
# MODEL IN STAN CODE
#-------------------------------------------------------------------------------
stan_code <- '
data {
  // Define variables in data
  // Number of level-1 observations (an integer)
  int<lower=0> Ni;
  // Number of level-2 clusters
  int<lower=0> Nj;
  // Number of level-3 clusters
  int<lower=0> Nk;

  // Cluster IDs
  int<lower=1> level2[Ni];
  int<lower=1> level3[Ni];
  int<lower=1> lev3ForLev2[Nj];

  // Continuous outcome
  real Y_ijk[Ni];
  // Continuous predictor
  real X_1ijk[Ni];
}

parameters {
  // Define parameters to estimate
  // Population intercept (a real number)
  real beta_0;
  // Population slope
  real beta_1;

  // Level-1
  real<lower=0> sigma_e0;

  // Level-2 random effect
  real u_0jk[Nj];
  real<lower=0> sigma_u0jk;

  // Level-3 random effect
  real u_0k[Nk];
  real<lower=0> sigma_u0k;

  // overdispersion parameter
  real<lower=0> tau;

}

transformed parameters  {
  // Varying intercepts
  real beta_0jk[Nj];
  real beta_0k[Nk];

  // Individual mean
  real mu[Ni];

  // Varying intercepts definition
  // Level-3 (level-3 random intercepts)
  for (k in 1:Nk) {
    beta_0k[k] <- beta_0 + u_0k[k];
  }
  // Level-2 (level-2 random intercepts)
  for (j in 1:Nj) {
    beta_0jk[j] <- beta_0k[lev3ForLev2[j]] + u_0jk[j];
  }
  // Individual mean
  // for (i in 1:Ni) {
    //mu[i] <- beta_0jk[level2[i]] + beta_1 * X_1ijk[i];
  //}
}

model {
  // Prior part of Bayesian inference
  // Flat prior for mu (no need to specify if non-informative)

  // Random effects distribution
  u_0k  ~ normal(0, sigma_u0k);
  u_0jk ~ normal(0, sigma_u0jk);
  tau   ~ gamma(0.001, 0.001);

  // Likelihood part of Bayesian inference
  // Outcome model N(mu, sigma^2) (use SD rather than Var)
  for (i in 1:Ni) {
  	
  	mu[i] ~ normal(beta_0jk[level2[i]] + beta_1 * X_1ijk[i], tau);
    Y_ijk[i] ~ poisson(exp(mu[i]));
      	
    // Y_ijk[i] ~ normal(mu[i], sigma_e0);
  }
}
'

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stan_results <- stan(model_code = stan_code, data = stan_data, chains = 8, iter = 10000, warmup = 500, thin = 10)

# note that the final number of iterations you'll be dealing with is:
# (( N_iterations - warmup ) / thin )*N_chains

get_sampler_params(stan_results)

stan_results_ext <- rstan::extract(stan_results, permuted = TRUE)

print(stan_results, pars = c("beta_0", "beta_1", "sigma_e0", "sigma_u0jk", "sigma_u0k"))

rstan::traceplot(stan_results, pars = c("beta_0","beta_1","sigma_e0","sigma_u0jk","sigma_u0k"), inc_warmup = FALSE)
beta_0
beta_1
sigma_e0
sigma_u0jk
sigma_u0k
str(stan_results_ext)
names(stan_results)
class(stan_results)


names(stan_results_ext)

class(stan_results_ext[["beta_0"]])

