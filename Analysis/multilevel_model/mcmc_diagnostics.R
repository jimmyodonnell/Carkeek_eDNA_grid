# important things:
# Obs v predicted 
# Acf 
# Prior against marginal posterior 
# Also traces 

# presumes your JAGS output is called "jags_out"

# and if you input:
# jags_params <- c(
                # "beta",
                # "sigma2",
                # "tau2",
                # "phi",
                # "P",
                # "beta_0"
# )

jags_params[6]
# check out the structure of the output
summary(jags_out)
str(jags_out)
names(jags_out) # "model" "BUGSoutput" "parameters.to.save" "model.file" "n.iter" "DIC"


jags_out$parameters.to.save[! jags_out$parameters.to.save %in% "deviance" ]

str(jags_out$BUGSoutput)
names(jags_out$BUGSoutput)


# simulations list
mcmc_list <- jags_out$BUGSoutput$sims.list
names(mcmc_list)
(mcmc_list[jags_params[6]])
mcmc_list[2]


# simulations array
mcmc_array <- jags_out$BUGSoutput$sims.array
str(mcmc_array)
dim(mcmc_array)
# mcmc_array[ iterations , chains , parameter ]
# N_iter
# N_chain

jags_params[]
names(mcmc_array[1,1,1])



param_current <- "sigma2"

xdim <- c(min(density(mcmc_array[,, param_current])$x), max(density(mcmc_array[,, param_current])$x))
ydim <- c(min(density(mcmc_array[,, param_current])$y), max(density(mcmc_array[,, param_current])$y))
plot(density(mcmc_array[,,param_current]), xlim = xdim, ylim = ydim, type = "l")
points(density(mcmc_array[,3,1]), type = "l")
xrange <- c(0,0.3)
plot("Density", xlim = xrange, ylim = ydim, type = "n", main = param_current)
points(density(mcmc_array[,,param_current]), type='l')
pred_x <- seq(xrange[1], xrange[2], 0.01)
points(pred_x, dgamma(pred_x, 0.01, 0.01), type='l', lty = 2)
legend("topright", legend = c("Marginal Posterior", "Prior"), lty = c(1,2), bty = "n")

min(mcmc_array[,, param_current])
boxplot(
	list(
		mcmc_array[,, param_current],
		dgamma(pred_x, 0.01, 0.01)
		)
	)


par(mfrow = c(4,1))
for(i in 1:dim(mcmc_array)[2]){
	acf(mcmc_array[,i, "sigma2"], main = paste("Estimate of sigma^2, chain", i))
}

rgamma()

xv = seq(0,10,.01)
plot(xv, dgamma(xv, 0.01, 0.01), type='l')

# simulations matrix: a matrix with (N_iter*N_chain) rows, and a column per parameter
dimnames(jags_out$BUGSoutput$sims.matrix)
class(jags_out$BUGSoutput$sims.matrix)
dim(jags_out$BUGSoutput$sims.matrix)
mcmc_matrix <- jags_out$BUGSoutput$sims.matrix
plot(density(mcmc_matrix[,15]))



jags_out$BUGSoutput$pD
jags_out$BUGSoutput$DIC

dim(jags_out$BUGSoutput$sims.array)
dimnames(jags_out$BUGSoutput$sims.array)
dimnames(jags_out$BUGSoutput$sims.array)[[3]] # parameter names

dim(jags_out$BUGSoutput$sims.array[,,])




mcmc_summary <- jags_out$BUGSoutput$summary



################################################################################################
# PLOT TRACES
################################################################################################

mcmc_array <- jags_out$BUGSoutput$sims.array

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
names(jags_out$BUGSoutput$sims.list)
# plot(jags_out$sims.list$population.mean, jags_out$sims.list$population.sd)

dim(jags_out$BUGSoutput$sims.list$deviance) # N_iter*N_chain by 1


boxplot(
	x = jags_out$BUGSoutput$sims.list$deviance,
	main = "Deviance"
	)

dim(jags_out$BUGSoutput$sims.list[[1]]) # matrix of dim 300000 (N_iter*N_chain) by 10 (N_taxa)

# so, to plot the estimated proportion of taxon 1 from all samples:
boxplot(jags_out$BUGSoutput$sims.list[[1]][,1])

# or, plot the estimates of P for all species:

pdf(file = file.path(fig_dir, "P_boxplot.pdf"))
boxplot(
	x = jags_out$BUGSoutput$sims.list[[1]],
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
