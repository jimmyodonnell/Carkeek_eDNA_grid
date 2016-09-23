
# INPUT COMES FROM SCRIPT "load_data.r"
# Watch out: if number of clusters is >9, the plotting will look funky
# Only single digits can be used as symbols in base plotting
#-------------------------------------------------------------------------------

EXPORT <- FALSE

library(fpc) #pamk()
library(cluster) #pam
library(vegan) # vegdist()

# choose between using proportional or raw counts
my_table <- otu_table[["mean"]] # otu_mean, otu_spvar, otu_named, as.binary(otu_mean), [,1:100] , otu_filt
my_metadata <- metadata[["mean"]] #metadata[!duplicated(metadata[,"env_sample_name"]),]

mydist <- vegdist(my_table, method = "bray", binary = FALSE) # , binary = TRUE

# or use `comm_dist` list of distance matrices (from distance_decay.R)

max_clusters <- nrow(my_table)-1

mypam <- pamk(data = mydist, 
              # criterion="ch", 
              krange = 1:(attributes(mydist)$Size-1)
              # krange = 6
              ) # to restrict range of Ks considered: , krange = 2:4

# run pam at each of a set of K (don't let it chose the best)
pam_list <- list()
for(i in 1:max_clusters){
  pam_list[[i]] <- pamk(data = mydist, krange = i)
}

# extract the silhouette data
sils <- lapply(
  pam_list[2:length(pam_list)], 
  function(x) silhouette(x$pamobject)
)
names(sils) <- as.character(2:length(pam_list))

# extract just the silhouette widths
sil_width <- lapply(
  sils, function(x) x[ , "sil_width"]
)

if(EXPORT){
  plot_base   <- "pam_sil"
  pdf_file    <- file.path(fig_dir, paste(plot_base, ".pdf", sep = ""))
  legend_file <- file.path(fig_dir, paste(plot_base, "_legend.txt", sep = ""))
  writeLines(
"Silhouette widths of PAM analysis. Points are the width of the PAM silhouette of each sample at each number of clusters (K). Red line is the mean, blue line is the median. Boxes encompass the interquartile range with a line at the median, and the whiskers extend to 1.5 times the interquartile range. Boxplot outliers are omitted for clarity.",
  con = legend_file)
  pdf(file = pdf_file, width = 8, height = 4) #, width = 8, height = 3
}
par(mar = c(4,4,1,1))
boxplot(sil_width, 
  outpch = NA, # suppress plotting of outliers
  xaxt = "n",
  las  = 1, 
  xlab = "Number of Clusters", ylab = "Silhouette Width")
xaxis_lab <- names(sil_width)
xaxis_lab[seq(from = 2, to = length(xaxis_lab), by = 2)] <- ""
axis(1, at = seq_along(sil_width), labels = xaxis_lab)
stripchart(sil_width, 
  method = "jitter", 
  pch = 21, col = "maroon", bg = hsv(1,0.5,1, alpha = 0.2), 
  add = TRUE, 
  vertical = TRUE)
lines(sapply(sil_width, mean), col = "red", lwd = 3)
lines(sapply(sil_width, median), col = "blue", lwd = 3)
if(EXPORT){
  dev.off()
}

# alternate plot:
# plot(sapply(sil_width, mean), type = "b")

pam_out <- mypam$pamobject
(K_optim <- mypam$nc) # number of clusters

mycolors <- gghue(mypam$nc)

dev.new()
# see ?clusplot.default
# pdf(file = file.path(fig_dir, "pam_plot.pdf"))
	par(mar = c(5,5,3,1))
	clusplot(
		x = as.matrix(mydist),
		clus = pam_out$clustering,
		sub = NULL,
		lines = 0,
		labels = 4,
		color = FALSE,
		col.p = mycolors[pam_out$clustering],
		col.txt = mycolors, #levels(as.factor(mycolors[pam_out$clustering])),
		col.clus = mycolors, #levels(as.factor(mycolors[pam_out$clustering])),
		main = NA
	)
# dev.off()

dev.new(width = 4, height = 7)

# USING MEAN DATA
# pdf(file = file.path(fig_dir, "pam_in_space.pdf"), width = 5, height = 7)
plot(
		x = my_metadata[,colname_xcoord],
		xaxt = "n",
		xlim = c(-500, 2500),
		y = log(my_metadata[,colname_ycoord] + 100),
		yaxt = "n",
		# log = "y",
		col = mycolors[mypam$pamobject$clustering[my_metadata[,colname_env_sample]]],
		# bg = mycolors[mypam$pamobject$clustering[my_metadata$sample_id]],
		pch = 19, #as.character(mypam$pamobject$clustering[my_metadata$env_sample_name]-1)
		cex = 4,
		# cex = 1.5,
		main = paste("membership to PAM clusters, K =", K_optim),
		bty = "n", 
		las = 1,
		xlab = "Position along shore (meters)",
		ylab = "Position from 0 (meters)"
	)
text(
		x = my_metadata[,colname_xcoord],
		xlim = c(-500, 2500),
		y = log(my_metadata[,colname_ycoord] + 100),
		# log = "y",
		col = "white",
		# bg = mycolors[mypam$pamobject$clustering[my_metadata$sample_id]],
		labels = as.character(mypam$pamobject$clustering[my_metadata[,colname_env_sample]]),
		cex = 1.5,
		las = 1
	)
	axis(side = 2, 
		at = unique(log(my_metadata[,colname_ycoord] + 100)), 
		labels = unique(my_metadata[,colname_ycoord]), 
		line = -1, lwd = 0, lwd.ticks = 1, las = 1)
	axis(side = 1, 
		at = unique(my_metadata[,colname_xcoord]), 
		labels = unique(my_metadata[,colname_xcoord]), 
		line = 0, lwd = 0, lwd.ticks = 1, las = 1)
# dev.off()


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
