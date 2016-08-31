#!/usr/bin/env Rscript

# plot map of sampled site

tif_file <- ""

library(rgdal) # readGDAL()

tif_import <- readGDAL(tif_file)

class(tif_import)

summary(tif_import)

ignore_above_0 <- function(x){
# convert any values of a vector that are > 0 to 0
  x[ x > 0 ] <- 0
  return(x)
}

attributes(tif_import)$data$band1 <- ignore_above_0(attributes(tif_import)$data$band1)

attributes(tif_import)$data$band1 <- attributes(tif_import)$data$band1/10

out_dir   <- dirname(tif_file)
base_orig <- basename(tif_file)
out_base  <- gsub("\\.tif$", ".pdf", base_orig)

out_file <- file.path(out_dir, out_base)

pdf(file = out_file)
# par(mar = c(4,4,4,4))
plot(tif_import)
# axis(1)
dev.off()


