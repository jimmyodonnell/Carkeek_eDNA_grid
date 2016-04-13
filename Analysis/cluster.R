
# INPUT COMES FROM SCRIPT "load_data.r"

library(fpc) #pamk()
library(cluster) #pam
library(vegan) # vegdist()

# function for plot colors (sorta like ggplot)
gghue <- function(n){
	hues = seq(15, 375, length = n+1)
	hcl(h = hues, l = 65, c = 100)[1:n]
}

# choose between using proportional or raw counts
# DATA <- otu_table
# DATA <- otu_filt[,]
DATA <- otu_mean
metadata <- metadata_exp[!duplicated(metadata_exp[,"env_sample_name"]),]

mydist <- vegdist(DATA, method = "bray") # , binary = TRUE

mypam <- pamk(data = mydist) # to restrict range of Ks considered: , krange = 2:4
pam_out <- mypam$pamobject
mypam$nc # number of clusters

plot(pam(mydist, k = mypam$nc), which.plots = 1)

mycolors <- gghue(mypam$nc)

# see ?clusplot.default
pdf(file = file.path(fig_dir, "pam_plot.pdf"))
	par(mar = c(5,5,1,1))
	clusplot(
		x = as.matrix(mydist), 
		clus = pam_out$clustering, 
		sub = NULL, 
		lines = 0, 
		labels = 4, 
		color = FALSE, 
		col.p = mycolors[pam_out$clustering], 
		col.txt = mycolors, 
		col.clus = mycolors, 
		main = NA
	)
dev.off()


set.seed(2) # to control jitter
# pdf(file = file.path(fig_dir, "medoids_in_space.pdf"), width = 7, height = 7)
	plot(
		x = metadata$dist_along_shore, 
		xaxt = "n", 
		xlim = c(-500, 2500), 
		y = log(metadata$dist_from_shore + 10), 
		yaxt = "n",
		# log = "y",
		col = mycolors[mypam$pamobject$clustering[metadata$sample_id]], 
		# bg = mycolors[mypam$pamobject$clustering[metadata$sample_id]], 
		pch = as.character(mypam$pamobject$clustering[metadata$sample_id]-1), 
		cex = 1, 
		# cex = 1.5, 
		main = "membership to PAM classifications", 
		las = 1, 
		xlab = "Position along shore (meters)", 
		ylab = "Position from 0 (meters)"
	)
	axis(side = 2, at = unique(log(metadata$dist_from_shore + 10)), labels = unique(metadata$dist_from_shore), las = 1)
	axis(side = 1, at = unique(metadata$dist_along_shore), labels = unique(metadata$dist_along_shore), las = 1)
# dev.off()


# USING MEAN DATA
pdf(file = file.path(fig_dir, "pam_in_space.pdf"), width = 5, height = 7)
plot(
		x = metadata$dist_along_shore, 
		xaxt = "n", 
		xlim = c(-500, 2500), 
		y = log(metadata$dist_from_shore + 10), 
		yaxt = "n",
		# log = "y",
		col = mycolors[mypam$pamobject$clustering[metadata$env_sample_name]], 
		# bg = mycolors[mypam$pamobject$clustering[metadata$sample_id]], 
		pch = 19, #as.character(mypam$pamobject$clustering[metadata$env_sample_name]-1)
		cex = 4, 
		# cex = 1.5, 
		main = "membership to PAM classifications", 
		las = 1, 
		xlab = "Position along shore (meters)", 
		ylab = "Position from 0 (meters)"
	)
points(
		x = metadata$dist_along_shore, 
		xaxt = "n", 
		xlim = c(-500, 2500), 
		y = log(metadata$dist_from_shore + 10), 
		yaxt = "n",
		# log = "y",
		col = "white", 
		# bg = mycolors[mypam$pamobject$clustering[metadata$sample_id]], 
		pch = as.character(mypam$pamobject$clustering[metadata$env_sample_name]), 
		cex = 2, 
		# cex = 1.5, 
		main = "membership to PAM classifications", 
		las = 1, 
		xlab = "Position along shore (meters)", 
		ylab = "Position from 0 (meters)"
	)
	axis(side = 2, at = unique(log(metadata$dist_from_shore + 10)), labels = unique(metadata$dist_from_shore), las = 1)
	axis(side = 1, at = unique(metadata$dist_along_shore), labels = unique(metadata$dist_along_shore), las = 1)
dev.off()