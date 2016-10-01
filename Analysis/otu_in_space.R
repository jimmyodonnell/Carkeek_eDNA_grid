# INPUT COMES FROM SCRIPT "load_data.r"

EXPORT <- FALSE

# REQUIRES
# 1. abundance table of categories (columns; e.g. "otus") that were observed (rows; e.g. "samples")
# 2. metadata file, with columns of:
#    - "lon" and "lat" columns, or x and y coordinates
#    - data matching rownames of abundance/otu table

my_table <- otu_table[["taxon"]] # otu_named # otu_filt , otu_mean
if(sum(colSums(my_table) < 1) > 0){
	print("removing otus with total abundance < 1")
	my_table <- my_table[,colSums(my_table) >= 1]
}
my_metadata <- metadata[["mean"]]

LOG_TRANSFORM <- FALSE
if(LOG_TRANSFORM) {
  my_table <- log(my_table + 1)
}

if(!identical(rownames(my_table), my_metadata[, colname_env_sample])){
	warning("Hold your horses: the rownames of the otu table don't match up with the metadata. This will make bad things happen.")
}

# plot scaled proportional abundance in space (i.e. darkest color is wherever OTU's proportional abundance is at maximum)

for(i in c(1:10)){
	plot(
		my_metadata[,colname_lon],
		my_metadata[,colname_lat],
		bg = rgb(r = 1, g = 0, b = 0, alpha = scale01(my_table[,i])),
		pch = 21,
		main = colnames(my_table)[i],
		xlab = "Longitude",
		ylab = "Latitude"
	)
}

# plot only a single otu, from
colnames(my_table)
select_taxa <- c(
  "Embiotocidae", # surfperch
  "Cupolaconcha meroclista", # vermetid
  "Thysanoessa raschii", # krill
  "Mytilus", # mussel
  "Sessilia" # barnacle
)
# "Cupolaconcha meroclista"
some_more <- c(
  "Chthamalus", # barnacle
  "Modiolus modiolus", # mussel
  "Mytilus trossulus", # mussel
  "Ptychodera flava", # very abundant Acorn worm
  "Panopea generosa", # geoduck
  "Sus scrofa", # pig
  "Dillwynella sp. SW-2012-1", # tiny snail
  "Thysanoessa", # krill
  "Sebastidae", # rockfish
  "Elysia pusilla", # slug
  "Homo sapiens" # humans
)
single_otu <- "Cupolaconcha meroclista" # 


# Plot in Lat/Lon space
# pdf(file = file.path(fig_dir, "otu_in_space_mean.pdf"), width = 7, height = 7)
plot(
	my_metadata[,colname_lon],
	my_metadata[,colname_lat],
	bg = rgb(r = 1, g = 0, b = 0, alpha = my_table[, single_otu]/max(my_table[, single_otu])),
	pch = 21,
	cex = my_table[, single_otu]/max(my_table[, single_otu])*1.5,
	# cex = 1.5,
	main = single_otu,
	xlab = "Longitude",
	ylab = "Latitude"
)
# add 'x' to points where abundance was 0
points(
	my_metadata[,colname_lon],
	my_metadata[,colname_lat],
	pch = ifelse(my_table[, single_otu]/max(my_table[, single_otu]) == 0, 4, NA_integer_)
)
# dev.off()


#-------------------------------------------------------------------------------
# Plot in from/alongshore space
#-------------------------------------------------------------------------------

plot_name <- "otu_in_space_select"

y_for_plotting <- as.numeric(as.factor(my_metadata[ , colname_ycoord]))

if(!exists("legend_text")){legend_text <- list()}

legend_text[plot_name] <- {
"Distribution of eDNA from select taxa.
Circles are colored and scaled by the proportion of that taxon's maximum proportional abundance.
That is, the largest circle is the same size in each of the panels, and occurs where that taxon contributed the greatest proportional abundance of reads to that sample."
}

if(EXPORT){
  pdf_file    <- file.path(fig_dir, paste(plot_name, ".pdf", sep = ""))
  legend_file <- file.path(fig_dir, paste(plot_name, "_legend.txt", sep = ""))
  writeLines(legend_text[[plot_name]], con = legend_file)
  pdf(file = pdf_file, width = 10, height = 4)
}

par(mfrow = c(1,length(select_taxa)))

for(i in 1:length(select_taxa)){
	point_data  <- scale01(my_table[, select_taxa[i]])
	point_color <- rgb(r = 1, g = 0, b = 0, alpha = point_data)
	point_size  <- point_data*5
	plot(
		x = my_metadata[ , colname_xcoord],
		xaxt = "n",
		xlim = c(-500, 2500), # add a little space to the right and left
		y = y_for_plotting,
		yaxt = "n",
		# log = "y",
		pch = 21,
		bg = point_color,
		cex = point_size,
		main = select_taxa[i],
		las = 1,
		xlab = "Distance along shore (meters)",
		ylab = "Distance from 0 (meters)"
	)
	points(
		x = my_metadata[ , colname_xcoord],
		y = y_for_plotting,
		# log = "y",
		cex = 1,
		pch = ifelse(
			point_data == 0,
			4, NA_integer_)
	)
    
	axis(side = 2, at = unique(y_for_plotting), labels = unique(my_metadata[,colname_ycoord]), las = 1)
	axis(side = 1, at = unique(my_metadata[,colname_xcoord]), labels = unique(my_metadata[,colname_xcoord]), las = 1)
}
if(EXPORT){
  dev.off()
}


#-------------------------------------------------------------------------------
# LOOP TO PRINT MANY AT THE SAME TIME
#-------------------------------------------------------------------------------

# Plot in from/alongshore space

plot_name <- "otu_in_space_50_named"

y_for_plotting <- as.numeric(as.factor(my_metadata[ , colname_ycoord]))
# old: log(my_metadata[ , colname_ycoord] + 10)

# if(!exists("legend_text")){legend_text <- list()}
# legend_text[plot_name] <- {"legend text goes here"}

if(EXPORT){
  pdf_file    <- file.path(fig_dir, paste(plot_name, ".pdf", sep = ""))
  # legend_file <- file.path(fig_dir, paste(plot_name, "_legend.txt", sep = ""))
  # writeLines(legend_text[[plot_name]], con = legend_file)
  pdf(file = pdf_file, width = 3, height = 7) #, width = 8, height = 3
}
for(i in 1:50){
	plot(
		x = my_metadata[ , colname_xcoord],
		xaxt = "n",
		xlim = c(-500, 2500), # add a little space to the right and left
		y = y_for_plotting,
		yaxt = "n",
		# log = "y",
		pch = 21,
		bg = rgb(r = 1, g = 0, b = 0, alpha = scale01(my_table[, i])),
		cex = scale01(my_table[, i])*5,
		# cex = 1.5,
		main = colnames(my_table)[i],
		las = 1,
		xlab = "Distance along shore (meters)",
		ylab = "Distance from 0 (meters)"
	)
	points(
		x = my_metadata[ , colname_xcoord],
		y = y_for_plotting,
		# log = "y",
		cex = 1,
		pch = ifelse(
			my_table[, i] == 0,
			4, NA_integer_)
	)
    
	axis(side = 2, at = unique(y_for_plotting), labels = unique(my_metadata[,colname_ycoord]), las = 1)
	axis(side = 1, at = unique(my_metadata[,colname_xcoord]), labels = unique(my_metadata[,colname_xcoord]), las = 1)
}
if(EXPORT){
  dev.off()
}
