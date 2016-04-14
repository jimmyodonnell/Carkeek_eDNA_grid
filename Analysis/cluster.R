
# INPUT COMES FROM SCRIPT "load_data.r"
# Watch out: if number of clusters is >9, the plotting will look funky
# Only single digits can be used as symbols in base plotting

library(fpc) #pamk()
library(cluster) #pam
library(vegan) # vegdist()

# function for plot colors (sorta like ggplot)
gghue <- function(n){
	hues = seq(15, 375, length = n+1)
	hcl(h = hues, l = 65, c = 100)[1:n]
}

# choose between using proportional or raw counts
# my_table <- otu_table
# my_table <- otu_filt[,]
my_table <- otu_mean # otu_mean, otu_spvar
my_metadata <- metadata_mean #metadata_exp[!duplicated(metadata_exp[,"env_sample_name"]),]

mydist <- vegdist(my_table, method = "bray") # , binary = TRUE

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


# USING MEAN DATA
pdf(file = file.path(fig_dir, "pam_in_space.pdf"), width = 5, height = 7)
plot(
		x = my_metadata[,colname_xcoord],
		xaxt = "n",
		xlim = c(-500, 2500),
		y = log(my_metadata[,colname_ycoord] + 10),
		yaxt = "n",
		# log = "y",
		col = mycolors[mypam$pamobject$clustering[my_metadata[,colname_env_sample]]],
		# bg = mycolors[mypam$pamobject$clustering[my_metadata$sample_id]],
		pch = 19, #as.character(mypam$pamobject$clustering[my_metadata$env_sample_name]-1)
		cex = 4,
		# cex = 1.5,
		main = "membership to PAM classifications",
		las = 1,
		xlab = "Position along shore (meters)",
		ylab = "Position from 0 (meters)"
	)
points(
		x = my_metadata[,colname_xcoord],
		xaxt = "n",
		xlim = c(-500, 2500),
		y = log(my_metadata[,colname_ycoord] + 10),
		yaxt = "n",
		# log = "y",
		col = "white",
		# bg = mycolors[mypam$pamobject$clustering[my_metadata$sample_id]],
		pch = as.character(mypam$pamobject$clustering[my_metadata[,colname_env_sample]]),
		cex = 2,
		# cex = 1.5,
		main = "membership to PAM classifications",
		las = 1,
		xlab = "Position along shore (meters)",
		ylab = "Position from 0 (meters)"
	)
	axis(side = 2, at = unique(log(my_metadata[,colname_ycoord] + 10)), labels = unique(my_metadata[,colname_ycoord]), las = 1)
	axis(side = 1, at = unique(my_metadata[,colname_xcoord]), labels = unique(my_metadata[,colname_xcoord]), las = 1)
dev.off()


# PLOT INDIVIDUAL PCR REPLICATES
if(nrow(my_table) > 24) {
	set.seed(2) # to control jitter
# pdf(file = file.path(fig_dir, "medoids_in_space.pdf"), width = 7, height = 7)
	plot(
		x = my_metadata[,colname_xcoord],
		xaxt = "n",
		xlim = c(-500, 2500),
		y = log(my_metadata[,colname_ycoord] + 10),
		yaxt = "n",
		# log = "y",
		col = mycolors[mypam$pamobject$clustering[my_metadata[,colname_sampleid]]],
		# bg = mycolors[mypam$pamobject$clustering[my_metadata$sample_id]],
		pch = as.character(mypam$pamobject$clustering[my_metadata[,colname_sampleid]]-1),
		cex = 1,
		# cex = 1.5,
		main = "membership to PAM classifications",
		las = 1,
		xlab = "Position along shore (meters)",
		ylab = "Position from 0 (meters)"
	)
	axis(side = 2, at = unique(log(my_metadata[,colname_ycoord] + 10)), labels = unique(my_metadata[,colname_ycoord]), las = 1)
	axis(side = 1, at = unique(my_metadata[,colname_xcoord]), labels = unique(my_metadata[,colname_xcoord]), las = 1)
# dev.off()
}
