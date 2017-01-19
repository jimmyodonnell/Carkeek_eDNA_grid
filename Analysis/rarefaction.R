# Calculate alpha diversity (species richness) and plot species accumulation curve
library(vegan)

EXPORT <- FALSE

dat      <- otu_table[["clean"]]
metadat  <- metadata[["clean"]]

# sum the counts of each OTU across PCR replicates within a sample
dat_reduced <- do.call(rbind, 
  lapply(
    split(as.data.frame(dat), metadat[,colname_env_sample]), 
    colSums)
)

dat <- dat_reduced

# change rownames for plotting
rownames(dat) <- gsub("PCT-", "", rownames(dat))

min_reads <- min(rowSums(dat))

alpha_div <- rarefy(dat, sample = min_reads)

rarefied  <- rrarefy(dat, sample = min_reads)

plot_name <- "rarefaction"

if(!exists("legend_text")){legend_text <- list()}
legend_text[plot_name] <- {
"Accumulation of OTUs from 24 environmental samples using randomized rarefaction.
Four replicate PCRs were conducted using DNA each environmental sample and independently sequenced, but these are collapsed here to illustrate a single representation of richness.
Sample names indicate the position in the sampling grid: south (S), central (C), or north (N), followed by the distance along the transect, in meters (0, 75, 125, 250, 500, 1000, 2000, 4000).
Vertical line indicates the minimum combined number of sequence reads per sample.
Horizontal lines indicate OTU richness for each sample at the minimum combined number of sequence reads."
}

if(EXPORT){
  pdf_file    <- file.path(fig_dir, paste(plot_name, ".pdf", sep = ""))
  legend_file <- file.path(fig_dir, paste(plot_name, "_legend.txt", sep = ""))
  writeLines(legend_text[[plot_name]], con = legend_file)
  pdf(file = pdf_file, width = 5, height = 4) #, width = 8, height = 3
}

par(mar = c(4,4,1,1))
rarecurve(dat, step = 10000, sample = min_reads, las = 1, 
  xlab = "Number of sequences", ylab = "OTUs", col = "orchid")

if(EXPORT){
  dev.off()
}

