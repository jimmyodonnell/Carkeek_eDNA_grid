
# MCMC diagnostics
# source http://www.johnmyleswhite.com/notebook/2010/08/29/mcmc-diagnostics-in-r-with-the-coda-package/

library('rjags')
 
N <- 1000
x <- 1:N
epsilon <- rnorm(N, 0, 1)
y <- x + epsilon

my_jags_model <- "
model
{
	for (i in 1:N)
	{
		y[i] ~ dnorm(y.hat[i], tau)
		y.hat[i] <- a + b * x[i]
	}
	a ~ dnorm(0, .0001)
	b ~ dnorm(0, .0001)
	tau <- pow(sigma, -2)
	sigma ~ dunif(0, 100)
}
"

jags <- jags.model(file = textConnection(my_jags_model),
                   data = list('x' = x,
                               'y' = y,
                               'N' = N),
                   n.chains = 4,
                   n.adapt = 10)
                   
samples <- coda.samples(jags,
                        c('a', 'b'),
                        1000)
                        

plot(samples)




jags <- jags.model(file = textConnection(my_jags_model),
                   data = list('x' = x,
                               'y' = y,
                               'N' = N),
                   n.chains = 4,
                   n.adapt = 1000)


samples <- coda.samples(jags,
                        c('a', 'b'),
                        1000)
                        
plot(samples)

gelman.plot(samples)

# Unfortunately, given our current call to jags.model(), it’s quite hard visually to identify convergence using Gelman plots, since the scales of these plots are not identical across our two examples, and the most prominent visual patterns are likely to be the results of random noise. There is a reason for this difficulty: we’re not properly initializing our sampler’s starting values separately for each chain. Both chains start from identical positions, which means that we don’t have enough power to really see the size of the space a theoretical chain might pass through before settling down. To fix that, we change our call to jags.model() to include an inits value, for which we provide an anonymous function that provides random values consistent with the prior we specified in example.bug. First, let’s repeat our previous approach again:

jags <- jags.model(file = textConnection(my_jags_model),
                   data = list('x' = x,
                               'y' = y,
                               'N' = N),
                   inits = function ()
                   {
                     list('a' = rnorm(1, 0, 100),
                          'b' = rnorm(1, 0, 100))
                   },
                   n.chains = 4,
                   n.adapt = 1000)
                   
                   
samples <- coda.samples(jags,
                        c('a', 'b'),
                        1000)

plot(samples)

gelman.plot(samples)
