pdf(file = file.path(fig_dir, "diss_by_dist_point.pdf"))

for ( i in 1:24) {
	
	# i <- 9

	geo <- as.matrix(geo_dist)[i,]
	geo <- geo[ geo > 0]
	
	com <- as.matrix(comm_dist)[i,]
	com <- com[ com > 0]
	
	the_name <- rownames(as.matrix(comm_dist))[i]

	# fit michaelis menten curve
	the_data <- data.frame(comm = com, space = geo)
	mm_fit <- nls(
	  formula = comm ~ Vm * space/(Km + space),
	  data = the_data,
	  start = list(
	    Vm = max(the_data$comm), # comment out and set Vm to 1 if asymptote should be 1
	    Km = max(the_data$comm)/2
	  )
	)
	mm_prediction <- predict(mm_fit)
	
	geo_dist_scaled <- log(geo_dist + 100)
	plot_x <- geo # geo_dist_scaled
	
	# pdf(file = file.path(fig_dir, "diss_by_dist.pdf")) #, width = 8, height = 3
		par(mar = c(4,5,1,1))
		plot(
			x = c(0, plot_x),
			y = c(0, com),
			ylim = c(0,1),
			xlim = c(0,4000), 
			xaxt = "n",
			pch = 21,
			las = 1,
			cex = 1,
			col = 1,
			bg = rgb(0,0,0,alpha = 0.1 ), #,alpha = 0.1
			xlab = "Distance between samples (meters)",
			ylab = "Bray-Curtis dissimilarity", 
			# log = "x"
		)
		axis(side = 1)
		abline(h = mean(com))
		#, at = unique(log(metadata$dist_from_shore + 100)), labels = unique(metadata$dist_from_shore)
	# abline(v = unique(log(metadata$dist_from_shore + 100)))
	
	# Add Michaelis Menten Fit
	lines(x = c(0, sort(geo)), y = c(0, sort(mm_prediction)), col = "red", lwd = 2, lty = 2)
	
	legend("bottomright", legend = the_name, bty = "n")
}

dev.off()
