
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

