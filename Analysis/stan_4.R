data_dir <- file.path("..", "Data")

data_file <- file.path(data_dir, "data_full_long_20.csv")

data_full <- read.csv(data_file)

mydata <- data_full

colname_level1 <- "sequenced_sample"
colname_level2 <- "transect_line"
colname_level3 <- "transect_position"
colname_level4 <- "OTU"
colname_predictor <- "dist_from_shore"
colname_outcome <- "count"

lookup_3in2 <- unique(mydata[,c(colname_level2, colname_level3)])[,colname_level3]
lookup_4in3 <- as.character(unique(mydata[,c(colname_level3, colname_level4)])[,colname_level4])

# plot(mydata[,colname_predictor], mydata[,colname_outcome])

stan_data <- list(Ni          = length(unique(mydata[,colname_level1])), 
                  Nj          = length(unique(mydata[,colname_level2])), 
                  Nk          = length(unique(mydata[,colname_level3])), 
                  Nl          = length(unique(mydata[,colname_level4])), 
                  level2      = mydata[,colname_level2],
                  level3      = mydata[,colname_level3],
                  level4      = mydata[,colname_level4],
                  lev3ForLev2 = lookup_3in2,
                  lev4ForLev2 = lookup_4in3,
                  Y_ijk       = mydata[,colname_outcome],
                  X_1ijk      = mydata[,colname_predictor]
                  )

lm_out <- lm(mydata[,colname_outcome] ~ mydata[,colname_predictor] )
summary(lm_out)