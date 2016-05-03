# INPUT COMES FROM SCRIPT "load_data.r"

# REQUIRES
# 1. abundance table of categories (columns; e.g. "otus") that were observed (rows; e.g. "samples")
# 2. metadata file, with columns of:
#    - "lon" and "lat" columns, or x and y coordinates
#    - data matching rownames of abundance/otu table

my_table <- otu_filt # otu_named # otu_filt , otu_mean
if(sum(colSums(my_table) < 1) > 0){
	print("removing otus with total abundance < 1")
	my_table <- my_table[,colSums(my_table) >= 1]
}
my_metadata <- metadata_mean

LOG_TRANSFORM <- FALSE
if(LOG_TRANSFORM) {
  my_table <- log(my_table + 1)
}

if(!identical(rownames(my_table), my_metadata[, colname_env_sample])){
	warning("Hold your horses: the rownames of the otu table don't match up with the metadata. This will make bad things happen.")
}

# this function rescales a numeric vector to 0 and 1
scale01 <- function(x){(x-min(x))/(max(x)-min(x))}

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

# plot only a single otu
single_otu <- "DUP_3"
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

set.seed(8) # to control jitter

################################################################################
# LOOP TO PRINT MANY AT THE SAME TIME

pdf(file = file.path(fig_dir, "otu_in_space_50.pdf"), width = 7, height = 7)
for(i in 1:50){
	plot(
		x = my_metadata[ , colname_xcoord],
		xaxt = "n",
		xlim = c(-500, 2500),
		y = log(my_metadata[ , colname_ycoord] + 10),
		yaxt = "n",
		# log = "y",
		bg = rgb(r = 1, g = 0, b = 0, alpha = scale01(my_table[, i])),
		pch = 21,
		cex = scale01(my_table[, i])*5,
		# cex = 1.5,
		main = colnames(my_table)[i],
		las = 1,
		xlab = "Position along shore (meters)",
		ylab = "Position from 0 (meters)"
	)
	points(
		x = my_metadata[ , colname_xcoord],
		y = log(my_metadata[ , colname_ycoord] + 10),
		# log = "y",
		cex = 1,
		pch = ifelse(
			my_table[, i] == 0,
			4, NA_integer_)
	)
	axis(side = 2, at = unique(log(my_metadata[,colname_ycoord] + 10)), labels = unique(my_metadata[,colname_ycoord]), las = 1)
	axis(side = 1, at = unique(my_metadata[,colname_xcoord]), labels = unique(my_metadata[,colname_xcoord]), las = 1)
}
dev.off()
