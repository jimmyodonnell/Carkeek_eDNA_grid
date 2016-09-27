
# INPUT COMES FROM SCRIPT "load_data.r"
# Watch out: if number of clusters is >9, the plotting will look funky
# Only single digits can be used as symbols in base plotting
#-------------------------------------------------------------------------------

EXPORT <- FALSE

# initiate a list to store legend text if it doesn't already exist
if(!exists("legend_text")){
  legend_text <- list()
}

# library(fpc) #pamk()
library(cluster) #pam
library(vegan) # vegdist()

# choose between using proportional or raw counts
my_table    <- otu_table[["mean"]] # mean, spvar, mean_unfilt, scale01, [,1:100]
my_metadata <- metadata[["mean"]]

mydist <- vegdist(my_table, method = "bray", binary = FALSE) # , binary = TRUE

# or convert similarity to dissimilarity from `comm_dist`, 
# a list of distance matrices (from distance_decay.R)

## cluster using PAM (partitioning around medoids)
pam_max <- nrow(my_table)-1

## run pam at each of a set of K (don't let it chose the best)
pam_list <- list()
for(i in 1: pam_max){
  pam_list[[i]] <- pam(x = mydist, k = i)
}

## extract the silhouette widths
pam_sil <- lapply(
  pam_list[2:length(pam_list)], 
  function(x) silhouette(x)[, "sil_width"]
)
names(pam_sil) <- as.character(2:length(pam_list))
pam_choice <- which.max(sapply(pam_sil, mean))

## cluster using FANNY (fuzzy analysis clustering)
fanny_max <- nrow(my_table)/2-1
fanny_list <- list()
for(i in 1:fanny_max){
  fanny_list[[i]] <- fanny(x = mydist, k = i)
}

fanny_sil <- lapply(
  fanny_list[2:length(fanny_list)], # can't get a silhouette from K = 1
  function(x) silhouette(x)[, "sil_width"]
)
names(fanny_sil) <- as.character(2:length(fanny_list))

fanny_choice <- which.max(sapply(fanny_sil, mean))


# pick a silhouette output to plot
sil_width <- pam_sil # fanny_sil, pam_sil


plot_name   <- "pam_sil"
legend_text[plot_name] <- "Silhouette widths from PAM analysis. Points are the width of the PAM silhouette of each sample at each number of clusters (K). Red line is the mean, blue line is the median. Boxes encompass the interquartile range with a line at the median, and the whiskers extend to 1.5 times the interquartile range. Boxplot outliers are omitted for clarity."

if(EXPORT){
  pdf_file    <- file.path(fig_dir, paste(plot_name, ".pdf", sep = ""))
  legend_file <- file.path(fig_dir, paste(plot_name, "_legend.txt", sep = ""))
  writeLines(legend_text[[plot_name]], con = legend_file)
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

#-------------------------------------------------------------------------------
## Plot the clustering results in two dimensions, a la PCA/NMDS
## This seems to be generally uninformative
plot_name   <- "pam_2d"
legend_text[plot_name] <- "Two-dimensional representation of PAM results. This plot is not really informative."

if(EXPORT){
  pdf_file    <- file.path(fig_dir, paste(plot_name, ".pdf", sep = ""))
  legend_file <- file.path(fig_dir, paste(plot_name, "_legend.txt", sep = ""))
  writeLines(legend_text[[plot_name]], con = legend_file)
  pdf(file = pdf_file, width = 5, height = 5)
}

for(i in 1:8){
whichK <- i
mycolors <- gghue(whichK)
	par(mar = c(5,5,3,1))
	clusplot(
		x = as.matrix(mydist),
		clus = pam_list[[whichK]]$clustering,
		sub = NULL,
		lines = 0,
		labels = 4,
		color = FALSE,
		col.p = mycolors[pam_list[[whichK]]$clustering],
		col.txt = mycolors, #levels(as.factor(mycolors[pam_out$clustering])),
		col.clus = mycolors, #levels(as.factor(mycolors[pam_out$clustering])),
		main = NA
	)
	legend("topleft", bty = "n", title = paste("K =", i), cex = 2, legend = "")
}
if(EXPORT){
  dev.off()
}

#-------------------------------------------------------------------------------
## Plot the sites in space, colored by cluster, for ALL levels of K
plot_name <- "pam_in_space_all"
legend_text[plot_name] <- "Cluster membership of sampled sites. Distance from onshore starting point is log scaled. Sites are colored and labeled by their assignment to a cluster by PAM analysis."

if(EXPORT){
  pdf_file    <- file.path(fig_dir, paste(plot_name, ".pdf", sep = ""))
  legend_file <- file.path(fig_dir, paste(plot_name, "_legend.txt", sep = ""))
  writeLines(legend_text[[plot_name]], con = legend_file)
  pdf(file = pdf_file, width = 5, height = 7)
}

for(i in 1:8){

whichK <- i
cluster_vec <- pam_list[[whichK]]$clustering
mycolors <- gghue(whichK)
site_cols <- mycolors[cluster_vec[my_metadata[,colname_env_sample]]]
site_labs <- LETTERS[cluster_vec[my_metadata[,colname_env_sample]]]

plot(
	x    = my_metadata[,colname_xcoord],
	xaxt = "n",
	xlim = c(-500, 2500),
	y    = log(my_metadata[,colname_ycoord] + 100),
	yaxt = "n",
	# log  = "y",
	col  = site_cols,
	# bg = mycolors[mypam$pamobject$clustering[my_metadata$sample_id]],
	pch  = 19, #as.character(mypam$pamobject$clustering[my_metadata$env_sample_name]-1)
	cex  = 4,
	main = paste("membership to PAM clusters, K =", whichK),
	bty  = "n", 
	las  = 1,
	xlab = "Position along shore (meters)",
	ylab = "Position from 0 (meters)"
)
text(
	x    = my_metadata[,colname_xcoord],
	xlim = c(-500, 2500),
	y    = log(my_metadata[,colname_ycoord] + 100),
	# log  = "y",
	col  = "white",
	# bg   = mycolors[mypam$pamobject$clustering[my_metadata$sample_id]],
	labels = site_labs,
	cex  = 1.5,
	las  = 1
)
axis(side = 2, 
	at = unique(log(my_metadata[,colname_ycoord] + 100)), 
	labels = unique(my_metadata[,colname_ycoord]), 
	line = -1, lwd = 0, lwd.ticks = 1, las = 1)
axis(side = 1, 
	at = unique(my_metadata[,colname_xcoord]), 
	labels = unique(my_metadata[,colname_xcoord]), 
	line = 0, lwd = 0, lwd.ticks = 1, las = 1)
}
if(EXPORT){
  dev.off()
}

#-------------------------------------------------------------------------------
## Plot the sites in space, colored by cluster, for ALL levels of K
plot_name <- "pam_in_space"
legend_text[plot_name] <- "Cluster membership of sampled sites. Distance from onshore starting point is log scaled. Sites are colored and labeled by their assignment to a cluster by PAM analysis for number of clusters (K) chosen based on a priori expectations (2) and mean silhouette width (8)."

if(EXPORT){
  pdf_file    <- file.path(fig_dir, paste(plot_name, ".pdf", sep = ""))
  legend_file <- file.path(fig_dir, paste(plot_name, "_legend.txt", sep = ""))
  writeLines(legend_text[[plot_name]], con = legend_file)
  pdf(file = pdf_file, width = 10, height = 7)
}

K_plot <- c(2,8)

par(mfcol = c(1, length(K_plot)), mar = c(4,5,1,1))
for(i in K_plot){

whichK <- i
cluster_vec <- pam_list[[whichK]]$clustering
mycolors <- gghue(whichK)
site_cols <- mycolors[cluster_vec[my_metadata[,colname_env_sample]]]
site_labs <- LETTERS[cluster_vec[my_metadata[,colname_env_sample]]]

plot(
	x    = my_metadata[,colname_xcoord],
	xaxt = "n",
	xlim = c(-500, 2500),
	y    = log(my_metadata[,colname_ycoord] + 100),
	yaxt = "n",
	# log  = "y",
	col  = site_cols,
	# bg = mycolors[mypam$pamobject$clustering[my_metadata$sample_id]],
	pch  = 19, #as.character(mypam$pamobject$clustering[my_metadata$env_sample_name]-1)
	cex  = 4,
	main = "",
	bty  = "n", 
	las  = 1,
	xlab = "Position along shore (meters)",
	ylab = "Position from 0 (meters)"
)
text(
	x    = my_metadata[,colname_xcoord],
	xlim = c(-500, 2500),
	y    = log(my_metadata[,colname_ycoord] + 100),
	# log  = "y",
	col  = "white",
	# bg   = mycolors[mypam$pamobject$clustering[my_metadata$sample_id]],
	labels = site_labs,
	cex  = 1.5,
	las  = 1
)
axis(side = 2, 
	at = unique(log(my_metadata[,colname_ycoord] + 100)), 
	labels = unique(my_metadata[,colname_ycoord]), 
	line = -1, lwd = 0, lwd.ticks = 1, las = 1)
axis(side = 1, 
	at = unique(my_metadata[,colname_xcoord]), 
	labels = unique(my_metadata[,colname_xcoord]), 
	line = 0, lwd = 0, lwd.ticks = 1, las = 1)
}
if(EXPORT){
  dev.off()
}

#-------------------------------------------------------------------------------
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
