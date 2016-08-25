rm(list=ls())
### GO GET THE NECESSARY LIBRARIES

##### TILAPIA IS OTU 4 ####

library(gtools)
library(R2jags)

# Set this to the subdirectory Analysis within the project folder (which should also contain subdirectories called Data and Figures)
# analysis_dir <- "/Users/threeprime/Documents/GoogleDrive/Kelly_Lab/Projects/Carkeek_eDNA_grid/Analysis"


# data_dir <- file.path("..", "Data")
# fig_dir <- file.path("..", "Figures")

setwd("/Users/threeprime/Documents/GoogleDrive/Kelly_Lab/Projects/Carkeek_eDNA_grid/Analysis/general_framework/Data")







#------------------------------------------
## Read in the clusters data + sort to make long form
# this file is a matrix of counts (integers) where
# rows = taxa
# columns = samples
# AND
# columns = distance measure from taxon's primer binding site and each of the tags (primer indexes)
# for original modeling paper, dim =  9,151; only top 10 OTUs
dat		<-	read.csv("abundance_distance_matrix+OLE.csv")
dat		<-	dat[dat$X!='OTU_4',] # remove Tilapia # no need once tables are separated appropriately

# OTU_table.csv
# primer_mismatch_both.csv
# primer_mismatch_F.csv
# primer_mismatch_R.csv
#------------------------------------------









#------------------------------------------
# Read in the "duplicate table":
# this file is a matrix of counts (integers)
# rows = unique DNA sequences ("DUP")
# columns = samples
# dim = (27973,75); 75 samples, 27973 unique sequences
dat.total.DNA	<-	read.csv("all_clusters.csv", header=TRUE, row.names = 1)

# duplicate_table.csv
#------------------------------------------

# calculate the total number of sequence counts per sample
dat.DNA			<-	colSums(dat.total.DNA)







#------------------------------------------
# Read in the sequencing metadata
# contains primer index data, environmental sample name, etc
dat.surv	<-	read.csv("sample_data.csv")
# sequencing_metadata.csv
#------------------------------------------


# isolate the primer index sequence ("tag") and the sample name
dat.surv.trim		<-	dat.surv[,c('Tag_Sequence','Sample')]

# exclude the last few samples (in this case, anything beyond row 22 is a positive or negative control sample)
dat.surv.trim		<-	dat.surv.trim[1:22,]

# add a column for the time the sample was taken (in this case embedded in the sample name)
dat.surv.trim$Time	<-	substr(dat.surv.trim$Sample,8,13)








################################################################################

# isolate the names of the OTUs
OTU		<-	as.character(dat$X)

# isolate the primer match data (genetic distance between primer sequences and template for each OTU sequence)
PRIMER	<-	dat[,77:ncol(dat)]
rownames(PRIMER) <- OTU

# how many OTUs are there?
N.sp	<-	length(OTU)

# make a vector of all of the primer index (tag) sequences
tags	<-	unique(substr(colnames(dat[,2:76]),1,6))


# isolate the sequence counts from each of the top 10 OTUs in each of the samples
COUNTS	<-  dat[,2:76]

# transpose the counts and make a dataframe
COUNTS			<- data.frame(t(COUNTS))

# name the columns the OTU names
colnames(COUNTS) <- OTU

# calculate the total number of sequence counts of the OTUs per sample
tot.obs			<- rowSums(COUNTS)

# add a column containing the counts of DNA sequences in each sample that do not belong to one of the OTUs
COUNTS$Other	<-	dat.DNA[match(rownames(COUNTS),names(dat.DNA))]
COUNTS$Other	<- COUNTS$Other - tot.obs

### STANDARDIZE THE NUMBER OF OTUs read to 100,000
# COUNTS_ORIG <- COUNTS
# adjust <- rowSums(COUNTS[,1:length(OTU)]) / 100000
# COUNTS <- round(COUNTS / adjust)

# ! COUNTS is now a dataframe where rows are samples, and


# add a column of the tag sequence
COUNTS$Tag	<-	substr(rownames(COUNTS),1,6)

# add a column of time the sample was taken
COUNTS$Time	<-	dat.surv.trim$Time[match(COUNTS$Tag,as.character(dat.surv.trim$Tag_Sequence))]

# exclude rows whose collection time is NA
COUNTS		<-	COUNTS[is.na(COUNTS$Time)==F,]

# order them by collection time
COUNTS		<-	COUNTS[order(COUNTS$Time),]

# add a column of zeros
COUNTS$t.numb		<-	0










# So, moving into the modeling aspect, the pieces of data we need are:
# 1. PRIMER = a matrix of distances between primer and template sequence for each OTU (rownames)
# 2. COUNTS = a data frame, rows = samples, cols = (counts of OTU1:OTU9+others, tag sequence, time sampled, time sampled as factor)
# 3.
# 4.
# 5.
# 6.














################################################################################
##### Manipulate COVARIATE DATA

# Make a matrix containing columns of OTU name, index sequence (tag), direction, genetic distance, and centralized distance
dat.covar	<- NULL
for(i in 1:ncol(PRIMER)){
	temp	<-	cbind(
				OTU=rownames(PRIMER),
				tag=substr(colnames(PRIMER)[i],5,10),
				direction=substr(colnames(PRIMER)[i],12,15),
				Dist=PRIMER[,i]
				)

	dat.covar	<-	rbind(dat.covar,temp)
}

# convert from matrix to data frame
dat.covar	<-	data.frame(dat.covar)

# make the distance column numeric
dat.covar$Dist	<-	as.numeric(as.character(dat.covar$Dist))

# split into separate data frames for each direction
dat.F	<-	dat.covar[dat.covar$direction=="F",]
dat.R	<-	dat.covar[dat.covar$direction=="R",]
dat.both	<-	dat.covar[dat.covar$direction=="both",]

# Add a column of centralized distance by subtracting the mean
dat.F$cent.dist		<-	dat.F$Dist - mean(dat.F$Dist)
dat.R$cent.dist		<-	dat.R$Dist - mean(dat.R$Dist)
dat.both$cent.dist	<-	dat.both$Dist - mean(dat.both$Dist)

# combine into an additional data frame
dat.COVAR	<-	data.frame(rbind(dat.F,dat.R,dat.both))




















#########################################################################################
#### ----- Define the covariates to use ------------------------------
#########################################################################################

# covariate_colname <- "direction"
# COVAR <- unique(dat.COVAR[,"direction"]) # ? as.character()
COVAR 	<-	c("F","R","both")






##### Turn Data and Covariates into arrays that can be used easily by JAGS

# the number of sampling events (11)
N.t	<- length(unique(COUNTS$Time))
# N.t		<-	N_times

for(i in 1:N.t){
	COUNTS$t.numb[COUNTS$Time == sort(unique(COUNTS$Time))[i]]	<-	i
}

colnames(COUNTS)[1:length(OTU)]	<-	OTU

# JO: I think the above code is another way of doing this:
t_numb <- as.numeric(as.factor(COUNTS$Time))
identical(t_numb, COUNTS$t.numb) # it is (as far as I understand the intent of the code)





# calculate the number of tags (for sampling event 1... assumed to be the same for all)
N.tag	<-	length(unique(COUNTS$Tag[COUNTS$t.numb==1]))

# more robust:
# N.tag <- unique(
	# sapply(
		# X = split(COUNTS$Tag, COUNTS$t.numb),
		# FUN = function(x) length(levels(as.factor(x)))
		# )
	# )
# if(length(N.tag) > 1)
	# stop("Different numbers of tags per sample!")



# calculate the number of replicates per tag
N.rep	<-	length(COUNTS$Tag[COUNTS$t.numb==1 & COUNTS$Tag == unique(COUNTS$Tag)[1]])

# more robust: N.rep <- unique(table(COUNTS$Tag))
# if(length(N.rep) > 1)
	# stop("Different numbers of reps per tag!")





# Make four dimensional arrays 	( # sequencer replicates, # OTUs, # Tags,  # Sample Times)
# JO: What is this array for?
DATA		<- 	array(
				data = 0,
				dim = c(N.rep, N.sp, N.tag, N.t),
				dimnames = list(
								REP = c("X1","X2","X3"),
								OTU = OTU,
								Tag = 1:N.tag,
								N.t = 1:N.t)
				)

# JO: What is this array for?
LABS		<-	array(
				data = 0,
				dim = c(N.rep,1,N.tag, N.t),
				dimnames = list(
								REP = c("X1","X2","X3"),
								NULL,
								Tag = 1:N.tag,
								N.t = 1:N.t
								)
				)

# JO: What is this array for?
# At this point, this array is identical to the array "DATA"
COVAR.mat	<-	array(
					data = 0,
					dim = c(N.rep,N.sp,N.tag, N.t),
					dimnames = list(
									REP = c("X1","X2","X3"),
									OTU = OTU,
									Tag = 1:N.tag,
									N.t = 1:N.t
								)
				)



# Fill in the arrays
for(i in 1:N.t){
		temp.tag 	<- unique(COUNTS$Tag[COUNTS$t.numb==i])
	for(j in 1:N.tag){
		temp		<- as.matrix(COUNTS[COUNTS$t.numb==i & COUNTS$Tag == temp.tag[j],1:N.sp])
		DATA[,,j,i]	<- temp
		LABS[,,j,i]	<- temp.tag[j]
	}
}

DATA[ 1, 9, 2, 11]

### THIS MAKES AN ARRAY FOR EACH COVARIATE BASE ON THE OBJECT "COVAR".
###  The arrays are named "COVAR.(VARIABLE NAME)"

for(m in 1:length(COVAR)){
	name.temp	<-	paste("COVAR",COVAR[m],sep=".")
	temp		<- dat.COVAR[dat.COVAR$direction == COVAR[m],]

	covar		<-	COVAR.mat
	for(i in 1:N.rep){
		for(j in 1:N.sp){
			for(k in 1:N.tag){
				for(l in 1:N.t){
					covar[i,j,k,l]	<-	temp$cent.dist[temp$OTU == dimnames(covar)$OTU[j]
											& temp$tag == LABS[i,1,k,l]]
				}
			}
		}
	}
	assign(name.temp,covar)
}


























##
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
#### ------------ START ESTIMATION ---------------------------------------
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

# In the following estimations, let:
# i = replicate (sequencing replicate in this case; 1:N.rep) (in the paper: "m")
# j = taxon (1:N.sp) (in the paper: "i")
# k = primer index per sample (1:N.tag)
# l = sampling event ("environmental sample" in time, in this case; 1:N.t)

# thus:
# DATA[i,j,k,l] = Z =  counts of sequences from replicate i for taxon j using primer index k from environmental sample l


# lambda =
# lambda[i,j,k,l] =

# beta = taxon-specific fixed effects ?( = fixed )?
# betas[j,l] =

# fixed =


# gammas = estimated coefficients (of covariation?)
# N.gammas = number of gamma parameters to use (I think this shouldn't be "user-defined" because it just depends on the number of covariates)
# thus, this is the parameter that defines how covariates shared across taxa (e.g. the quality of match between the primer and taxa DNA) will affect the observed number of DNA sequences for each taxon.

# COVAR.X = array for covariate "X"... known as "H" in the paper


# eta =
# eta[i,j,k,l] =

# sigma2 = stochasticity in the PCR and sequencing process

# phi =

# p =

# P =

# tau2 =
# varies uniformly among taxa















# Z = DATA
# H = COVAR.X

# general structure:
Z[i,j,k,l] ~ dpois(
				exp(
					betas[j,l] + gamma*H[i,j,k,l] + eta[i,j,k,l]
				)
			)
)










# 111111111111111111111111111111111111111111111111111111111111111111111111111111
##########################################################################################
#### BASE
##### FIXED for the match between the tags + primer and the substrate
#####	+ SIMPLEST RANDOM EFFECT
##########################################################################################
##########################################################################################


N.gammas	<-	1 #length(COVAR)

jagsscript = cat(
# "
model {
	for(l in 1:N.t){
	  	for(k in 1:N.tag){
		  	for(j in 1:N.sp){
	  			for(i in 1:N.rep){


  					DATA[i,j,k,l] ~ dpois(
  										exp(
  											lambda[i,j,k,l]
  										)
  									) # lambda is a vector of (non-negative integer) quantiles that define the Poisson distribution


				}
			}
 		}
 	}

	# Why do separate lambdas for the first taxon and everyone else?
	for(l in 1:N.t){
		for(k in 1:N.tag){
				for(i in 1:N.rep){

					fixed[i,1,k,l]		<- betas[1,l] 				+ gammas[1]*COVAR.F[i,1,k,l]# + gammas[2]*COVAR.R[i,1,k,l]
					lambda[i,1,k,l]		<- fixed[i,1,k,l] + eta[i,1,k,l]

				}
			# include taxon-specific effects for the other species?
  			for(j in 2:N.sp){
				for(i in 1:N.rep){
					fixed[i,j,k,l]		<- betas[1,l] + betas[j,l] 	+ gammas[1]*COVAR.F[i,j,k,l]# + gammas[2]*COVAR.R[i,j,k,l]
					lambda[i,j,k,l]		<- fixed[i,j,k,l] + eta[i,j,k,l]
				}
			}
		}
	}

	for(l in 1:N.t){
		for(k in 1:N.tag){
	  		for(j in 1:N.sp){
	  			for(i in 1:N.rep){
  					eta[i,j,k,l] ~ dnorm(0, 1/sigma2)
				}
 			}
 		}
	}

	### Derived Quantities
	for(l in 1:N.t){
		p[1,l] <- exp(betas[1,l])
		for(j in 2:N.sp){
			p[j,l]	<-	exp(betas[1,l] + betas[j,l])
		}
	}
	for(l in 1:N.t){
		for(j in 1:N.sp){
			P[j,l]	<-	p[j,l] / sum(p[,l])
		}
  	}

  	### Priors
	for(l in 1:N.t){
	  	for(j in 1:N.sp){
		 	betas[j,l] ~ dnorm(0,0.001)
    	}
	}
   	for(i in 1:N.gammas){
		 gammas ~ dnorm(0,0.001)
   	}

 	sigma2 ~ dunif(0,1000)

#  	for(j in 1:N.sp){
# 	 	tau2[j] ~ dunif(0,100)
# 	}


}
# "
, file="jags_dummy1.txt")

jags.data     = list("DATA","N.t","N.sp","N.tag","N.rep","COVAR.F","N.gammas")
jags.params   = c("betas","gammas","sigma2","P","fixed") #"tau2"
model.loc		= c("jags_dummy1.txt")

Nburn = 100000
Niter = 10000
Nchain=	3
Inits = NULL
for(i in 1:Nchain){
	Inits[[i]]	<- list("betas"=matrix(runif(N.t*N.sp,-10,10),N.sp,N.t),
 						"gammas"=runif(N.gammas,-10,10),
						"sigma2"=runif(1,1,10))
# 						"tau2"=runif(N.sp,1,10))
}
jags.model.rand.base = jags(
						data = jags.data,
						inits = Inits,
						parameters.to.save= jags.params,
						model.file=model.loc,
						n.chains = Nchain,
						n.burnin = Nburn,
						n.thin = 1,
						n.iter = Nburn+Niter,
						DIC = TRUE
)

















# 222222222222222222222222222222222222222222222222222222222222222222222222222222
##########################################################################################
## BASE 2Fix [ uses both the forward and the reverse match % as covariates ]

N.gammas	<-	2#length(COVAR)

jagsscript = cat(
#"
model {
	for(l in 1:N.t){
	  	for(k in 1:N.tag){
		  	for(j in 1:N.sp){
	  			for(i in 1:N.rep){
  					DATA[i,j,k,l] ~ dpois(exp(lambda[i,j,k,l]))
				}
			}
 		}
 	}

	for(l in 1:N.t){
		for(k in 1:N.tag){
				for(i in 1:N.rep){
					lambda[i,1,k,l]		<- betas[1,l] 				+ gammas[1]*COVAR.F[i,1,k,l] + gammas[2]*COVAR.R[i,1,k,l] + eta[i,1,k,l]
				}
  			for(j in 2:N.sp){
				for(i in 1:N.rep){
					lambda[i,j,k,l]		<- betas[1,l] + betas[j,l] 	+ gammas[1]*COVAR.F[i,j,k,l] + gammas[2]*COVAR.R[i,j,k,l] + eta[i,j,k,l]
				}
			}
		}
	}

	for(l in 1:N.t){
		for(k in 1:N.tag){
	  		for(j in 1:N.sp){
	  			for(i in 1:N.rep){
  					eta[i,j,k,l] ~ dnorm(0, 1/sigma2)
				}
 			}
 		}
	}

	### Derived Quantities
	for(l in 1:N.t){
		p[1,l] <- exp(betas[1,l])
		for(j in 2:N.sp){
			p[j,l]	<-	exp(betas[1,l] + betas[j,l])
		}
	}
	for(l in 1:N.t){
		for(j in 1:N.sp){
			P[j,l]	<-	p[j,l] / sum(p[,l])
		}
  	}

  	### Priors
	for(l in 1:N.t){
	  	for(j in 1:N.sp){
		 	betas[j,l] ~ dnorm(0,0.001)
    	}
	}
   	for(m in 1:N.gammas){
		 gammas[m] ~ dnorm(0,0.001)
    }


 	sigma2 ~ dunif(0,1000)

#  	for(j in 1:N.sp){
# 	 	tau2[j] ~ dunif(0,100)
# 	}


}
#",
file="jags_dummy1.txt")

jags.data     = list("DATA","N.t","N.sp","N.tag","N.rep","COVAR.F","COVAR.R","N.gammas") #
jags.params   = c("betas","gammas","sigma2","P") #"tau2"
model.loc		= c("jags_dummy1.txt")

Nburn = 100000
Niter = 10000
Nchain=	3
Inits = NULL
for(i in 1:Nchain){
	Inits[[i]]	<- list("betas"=matrix(runif(N.t*N.sp,-10,10),N.sp,N.t),
						"gammas"=runif(N.gammas,-10,10),
						"sigma2"=runif(1,1,10))
# 						"tau2"=runif(N.sp,1,10))
}
jags.model.rand.base.2fix = jags(jags.data, inits = Inits, parameters.to.save= jags.params, model.file=model.loc, n.chains = Nchain, n.burnin = Nburn, n.thin = 1, n.iter = Nburn+Niter, DIC = TRUE)















# 333333333333333333333333333333333333333333333333333333333333333333333333333333
##########################################################################################
## BASE Both

N.gammas	<-	1#length(COVAR)

jagsscript = cat(
#"
model {
	for(l in 1:N.t){
	  	for(k in 1:N.tag){
		  	for(j in 1:N.sp){
	  			for(i in 1:N.rep){
  					DATA[i,j,k,l] ~ dpois(exp(lambda[i,j,k,l]))
				}
			}
 		}
 	}

	for(l in 1:N.t){
		for(k in 1:N.tag){
				for(i in 1:N.rep){
					lambda[i,1,k,l]		<- betas[1,l] 				+ gammas[1]*COVAR.both[i,1,k,l]  + eta[i,1,k,l] #+ gammas[2]*COVAR.R[i,1,k,l]
				}
  			for(j in 2:N.sp){
				for(i in 1:N.rep){
					lambda[i,j,k,l]		<- betas[1,l] + betas[j,l] 	+ gammas[1]*COVAR.both[i,j,k,l]  + eta[i,j,k,l] #+ gammas[2]*COVAR.R[i,j,k,l]
				}
			}
		}
	}

	for(l in 1:N.t){
		for(k in 1:N.tag){
	  		for(j in 1:N.sp){
	  			for(i in 1:N.rep){
  					eta[i,j,k,l] ~ dnorm(0, 1/sigma2)
				}
 			}
 		}
	}

	### Derived Quantities
	for(l in 1:N.t){
		p[1,l] <- exp(betas[1,l])
		for(j in 2:N.sp){
			p[j,l]	<-	exp(betas[1,l] + betas[j,l])
		}
	}
	for(l in 1:N.t){
		for(j in 1:N.sp){
			P[j,l]	<-	p[j,l] / sum(p[,l])
		}
  	}

  	### Priors
	for(l in 1:N.t){
	  	for(j in 1:N.sp){
		 	betas[j,l] ~ dnorm(0,0.001)
    	}
	}
   	for(m in 1:N.gammas){
		 gammas[m] ~ dnorm(0,0.001)
    }


 	sigma2 ~ dunif(0,100)

#  	for(j in 1:N.sp){
# 	 	tau2[j] ~ dunif(0,100)
# 	}


}
#",
file="jags_dummy1.txt")

  jags.data     = list("DATA","N.t","N.sp","N.tag","N.rep","COVAR.both","N.gammas") #
  jags.params   = c("betas","gammas","sigma2","P") #"tau2"
  model.loc		= c("jags_dummy1.txt")

	Nburn = 100000
	Niter = 10000
	Nchain=	3
	Inits = NULL
	for(i in 1:Nchain){
		Inits[[i]]	<- list("betas"=matrix(runif(N.t*N.sp,-10,10),N.sp,N.t),
 						"gammas"=runif(N.gammas,-10,10),
						"sigma2"=runif(1,1,10))
# 						"tau2"=runif(N.sp,1,10))
	}
jags.model.rand.base.both = jags(jags.data, inits = Inits, parameters.to.save= jags.params, model.file=model.loc, n.chains = Nchain, n.burnin = Nburn, n.thin = 1, n.iter = Nburn+Niter, DIC = TRUE)





















# 444444444444444444444444444444444444444444444444444444444444444444444444444444
##########################################################################################
##########################################################################################
## BASE.C
##### FIXED for the match between the tags + primer and the substrate
#####	+ SINGLE RANDOM EFFECT one for all species but for each species-tag-time interaction (not each observation)
##########################################################################################

N.gammas	<-	1#length(COVAR)

jagsscript = cat(
#"
model {
	for(l in 1:N.t){
	  	for(k in 1:N.tag){
		  	for(j in 1:N.sp){
	  			for(i in 1:N.rep){
  					DATA[i,j,k,l] ~ dpois(exp(lambda[i,j,k,l]))
				}
			}
 		}
 	}

	for(l in 1:N.t){
		for(k in 1:N.tag){
				for(i in 1:N.rep){
					lambda[i,1,k,l]		<- betas[1,l] 				+ gammas*COVAR.F[i,1,k,l] + phi[1,k,l]
				}
  			for(j in 2:N.sp){
				for(i in 1:N.rep){
					lambda[i,j,k,l]		<- betas[1,l] + betas[j,l] 	+ gammas*COVAR.F[i,j,k,l] + phi[j,k,l]
				}
			}
		}
	}

# 	for(l in 1:N.t){
# 		for(k in 1:N.tag){
# 	  		for(j in 1:N.sp){
# 	  			for(i in 1:N.rep){
#   					eta[i,j,k,l] ~ dnorm(0, 1/sigma2)
# 				}
#  			}
#  		}
# 	}
	for(l in 1:N.t){
 		for(k in 1:N.tag){
	 		for(j in 1:N.sp){
 					phi[j,k,l]	~ dnorm(0, 1/tau2)
			}
		}
	}

	### Derived Quantities
	for(l in 1:N.t){
		p[1,l] <- exp(betas[1,l])
		for(j in 2:N.sp){
			p[j,l]	<-	exp(betas[1,l] + betas[j,l])
		}
	}
	for(l in 1:N.t){
		for(j in 1:N.sp){
			P[j,l]	<-	p[j,l] / sum(p[,l])
		}
  	}

  	### Priors
	for(l in 1:N.t){
	  	for(j in 1:N.sp){
		 	betas[j,l] ~ dnorm(0,0.001)
    	}
	}
#   	for(i in 1:N.gammas){
		 gammas ~ dnorm(0,0.001)
#    	}


#  	sigma2 	~ dunif(0,100)
 	tau2	~ dunif(0,1000)

}
#",
file="jags_dummy1.txt")

  jags.data     = list("DATA","N.t","N.sp","N.tag","N.rep","COVAR.F") #,"COVAR.R","N.gammas"
  jags.params   = c("betas","gammas","tau2","P") #,"sigma2"
  model.loc		= c("jags_dummy1.txt")

	Nburn = 100000
	Niter = 10000
	Nchain=	3
	Inits = NULL
	for(i in 1:Nchain){
		Inits[[i]]	<- list("betas"=matrix(runif(N.t*N.sp,-10,10),N.sp,N.t),
 						"gammas"=runif(N.gammas,-10,10),
# 						"sigma2"=runif(1,0,30),
 						"tau2"=runif(1,0,30))
	}
jags.model.base.c.F = jags(jags.data, inits = Inits, parameters.to.save= jags.params, model.file=model.loc, n.chains = Nchain, n.burnin = Nburn, n.thin = 1, n.iter = Nburn+Niter, DIC = TRUE)



























# 555555555555555555555555555555555555555555555555555555555555555555555555555555
##########################################################################################

N.gammas	<-	2#length(COVAR)

jagsscript = cat(
# "
model {
	for(l in 1:N.t){
	  	for(k in 1:N.tag){
		  	for(j in 1:N.sp){
	  			for(i in 1:N.rep){
  					DATA[i,j,k,l] ~ dpois(exp(lambda[i,j,k,l]))
				}
			}
 		}
 	}

	for(l in 1:N.t){
		for(k in 1:N.tag){
				for(i in 1:N.rep){
					lambda[i,1,k,l]		<- betas[1,l] 				+ gammas[1]*COVAR.F[i,1,k,l] + gammas[2]*COVAR.R[i,1,k,l] + phi[1,k,l]
				}
  			for(j in 2:N.sp){
				for(i in 1:N.rep){
					lambda[i,j,k,l]		<- betas[1,l] + betas[j,l] 	+ gammas[1]*COVAR.F[i,j,k,l] + gammas[2]*COVAR.R[i,1,k,l] + phi[j,k,l]
				}
			}
		}
	}

# 	for(l in 1:N.t){
# 		for(k in 1:N.tag){
# 	  		for(j in 1:N.sp){
# 	  			for(i in 1:N.rep){
#   					eta[i,j,k,l] ~ dnorm(0, 1/sigma2)
# 				}
#  			}
#  		}
# 	}
	for(l in 1:N.t){
 		for(k in 1:N.tag){
	 		for(j in 1:N.sp){
 					phi[j,k,l]	~ dnorm(0, 1/tau2)
			}
		}
	}

	### Derived Quantities
	for(l in 1:N.t){
		p[1,l] <- exp(betas[1,l])
		for(j in 2:N.sp){
			p[j,l]	<-	exp(betas[1,l] + betas[j,l])
		}
	}
	for(l in 1:N.t){
		for(j in 1:N.sp){
			P[j,l]	<-	p[j,l] / sum(p[,l])
		}
  	}

  	### Priors
	for(l in 1:N.t){
	  	for(j in 1:N.sp){
		 	betas[j,l] ~ dnorm(0,0.001)
    	}
	}
   	for(i in 1:N.gammas){
		 gammas[i] ~ dnorm(0,0.001)
   	}


#  	sigma2 	~ dunif(0,100)
 	tau2	~ dunif(0,1000)

}
# ",
file="jags_dummy1.txt")

jags.data     = list("DATA","N.t","N.sp","N.tag","N.rep","COVAR.F","COVAR.R","N.gammas")
jags.params   = c("betas","gammas","tau2","P") #,"sigma2"
model.loc		= c("jags_dummy1.txt")

Nburn = 100000
Niter = 10000
Nchain=	3
Inits = NULL
for(i in 1:Nchain){
	Inits[[i]]	<- list("betas"=matrix(runif(N.t*N.sp,-10,10),N.sp,N.t),
						"gammas"=runif(N.gammas,-10,10),
# 						"sigma2"=runif(1,0,30),
 						"tau2"=runif(1,0,30))
}

jags.model.base.c.F.R = jags(jags.data, inits = Inits, parameters.to.save= jags.params, model.file=model.loc, n.chains = Nchain, n.burnin = Nburn, n.thin = 1, n.iter = Nburn+Niter, DIC = TRUE)
































# 666666666666666666666666666666666666666666666666666666666666666666666666666666
##########################################################################################
##########################################################################################
## V2
##### FIXED for the match between the tags + primer and the substrate
#####	+ TWO RANDOM EFFECT one shared by all species and one for each species (with realizations for each tag and time)
##########################################################################################

N.gammas	<-	1#length(COVAR)

jagsscript = cat(
# "
model {
	for(l in 1:N.t){
	  	for(k in 1:N.tag){
		  	for(j in 1:N.sp){
	  			for(i in 1:N.rep){
  					DATA[i,j,k,l] ~ dpois(exp(lambda[i,j,k,l]))
				}
			}
 		}
 	}

	for(l in 1:N.t){
		for(k in 1:N.tag){
				for(i in 1:N.rep){
					lambda[i,1,k,l]		<- betas[1,l] 				+ gammas*COVAR.F[i,1,k,l] + eta[i,1,k,l] + phi[1,k,l] #+ gammas[2]*COVAR.R[i,1,k,l]
				}
  			for(j in 2:N.sp){
				for(i in 1:N.rep){
					lambda[i,j,k,l]		<- betas[1,l] + betas[j,l] 	+ gammas*COVAR.F[i,j,k,l] + eta[i,j,k,l] + phi[j,k,l] #+ gammas[2]*COVAR.R[i,j,k,l] + eta[i,j,k,l] + phi[i,j,k,l]
				}
			}
		}
	}

	for(l in 1:N.t){
		for(k in 1:N.tag){
	  		for(j in 1:N.sp){
	  			for(i in 1:N.rep){
  					eta[i,j,k,l] ~ dnorm(0, 1/sigma2)
				}
 			}
 		}
	}


	for(l in 1:N.t){
		for(k in 1:N.tag){
	 		for(j in 1:N.sp){
 					phi[j,k,l]	~ dnorm(0, 1/tau2[j])
			}
		}
	}



	### Derived Quantities
	for(l in 1:N.t){
		p[1,l] <- exp(betas[1,l])
		for(j in 2:N.sp){
			p[j,l]	<-	exp(betas[1,l] + betas[j,l])
		}
	}
	for(l in 1:N.t){
		for(j in 1:N.sp){
			P[j,l]	<-	p[j,l] / sum(p[,l])
		}
  	}

  	### Priors
	for(l in 1:N.t){
	  	for(j in 1:N.sp){
		 	betas[j,l] ~ dnorm(0,0.001)
    	}
	}
#   	for(i in 1:N.gammas){
		 gammas ~ dnorm(0,0.001)
#    	}


 	sigma2 ~ dunif(0,1000)

	  	for(j in 1:N.sp){
	 	 	tau2[j] ~ dunif(0,1000)
 		}


}
# ",
file="jags_dummy1.txt")

jags.data     = list("DATA","N.t","N.sp","N.tag","N.rep","COVAR.F","lambda") #,"COVAR.R","N.gammas"
jags.params   = c("betas","gammas","sigma2","tau2","P") #
model.loc		= c("jags_dummy1.txt")

Nburn = 100000
Niter = 10000
Nchain=	3
Inits = NULL
for(i in 1:Nchain){
	Inits[[i]]	<- list("betas"=matrix(runif(N.t*N.sp,-10,10),N.sp,N.t),
						"gammas"=runif(N.gammas,-10,10),
						"sigma2"=runif(1,0,30),
						"tau2"=runif(N.sp,0,30))
}
jags.model.rand.v2 = jags(jags.data, inits = Inits, parameters.to.save= jags.params, model.file=model.loc, n.chains = Nchain, n.burnin = Nburn, n.thin = 1, n.iter = Nburn+Niter, DIC = TRUE)
































# 777777777777777777777777777777777777777777777777777777777777777777777777777777
##########################################################################################
##########################################################################################
## V3
##### FIXED for the match between the tags + primer and the substrate
#####	+ TWO RANDOM EFFECT one shared by all species and one for each species interaction (with realizations for each time)
##########################################################################################

N.gammas	<-	1#length(COVAR)

jagsscript = cat(
# "
model {
	for(l in 1:N.t){
	  	for(k in 1:N.tag){
		  	for(j in 1:N.sp){
	  			for(i in 1:N.rep){
  					DATA[i,j,k,l] ~ dpois(exp(lambda[i,j,k,l]))
				}
			}
 		}
 	}

	for(l in 1:N.t){
		for(k in 1:N.tag){
				for(i in 1:N.rep){
					lambda[i,1,k,l]		<- betas[1,l] 				+ gammas*COVAR.F[i,1,k,l] + eta[i,1,k,l] + phi[1,l] #+ gammas[2]*COVAR.R[i,1,k,l]
				}
  			for(j in 2:N.sp){
				for(i in 1:N.rep){
					lambda[i,j,k,l]		<- betas[1,l] + betas[j,l] 	+ gammas*COVAR.F[i,j,k,l] + eta[i,j,k,l] + phi[j,l] #+ gammas[2]*COVAR.R[i,j,k,l] + eta[i,j,k,l] + phi[i,j,k,l]
				}
			}
		}
	}

	for(l in 1:N.t){
		for(k in 1:N.tag){
	  		for(j in 1:N.sp){
	  			for(i in 1:N.rep){
  					eta[i,j,k,l] ~ dnorm(0, 1/sigma2)
				}
 			}
 		}
	}


	for(l in 1:N.t){
	 		for(j in 1:N.sp){
 					phi[j,l]	~ dnorm(0, 1/tau2[j])
			}
	}



	### Derived Quantities
	for(l in 1:N.t){
		p[1,l] <- exp(betas[1,l])
		for(j in 2:N.sp){
			p[j,l]	<-	exp(betas[1,l] + betas[j,l])
		}
	}
	for(l in 1:N.t){
		for(j in 1:N.sp){
			P[j,l]	<-	p[j,l] / sum(p[,l])
		}
  	}

  	### Priors
	for(l in 1:N.t){
	  	for(j in 1:N.sp){
		 	betas[j,l] ~ dnorm(0,0.001)
    	}
	}
#   	for(i in 1:N.gammas){
		 gammas ~ dnorm(0,0.001)
#    	}


 	sigma2 ~ dunif(0,1000)
  	for(j in 1:N.sp){
 	 	tau2[j] ~ dunif(0,1000)
	}



}
# ",
file="jags_dummy1.txt")

jags.data     = list("DATA","N.t","N.sp","N.tag","N.rep","COVAR.F") #,"COVAR.R","N.gammas"
jags.params   = c("betas","gammas","sigma2","tau2","P") #
model.loc		= c("jags_dummy1.txt")

Nburn = 100000
Niter = 10000
Nchain=	3
Inits = NULL
for(i in 1:Nchain){
	Inits[[i]]	<- list("betas"=matrix(runif(N.t*N.sp,-10,10),N.sp,N.t),
						"gammas"=runif(N.gammas,-10,10),
						"sigma2"=runif(1,0,30),
						"tau2"=runif(N.sp,0,30))
}

jags.model.rand.v3 = jags(jags.data, inits = Inits, parameters.to.save= jags.params, model.file=model.loc, n.chains = Nchain, n.burnin = Nburn, n.thin = 1, n.iter = Nburn+Niter, DIC = TRUE)

































# 888888888888888888888888888888888888888888888888888888888888888888888888888888
##########################################################################################
##########################################################################################
## V3.B
##### FIXED for the match between the tags + primer and the substrate
#####	+ TWO RANDOM EFFECT one shared by all species (for each observation) and one for all species but for each species-tag-time interaction
##########################################################################################

N.gammas	<-	1#length(COVAR)

jagsscript = cat(
# "
model {
	for(l in 1:N.t){
	  	for(k in 1:N.tag){
		  	for(j in 1:N.sp){
	  			for(i in 1:N.rep){
  					DATA[i,j,k,l] ~ dpois(exp(lambda[i,j,k,l]))
				}
			}
 		}
 	}

	for(l in 1:N.t){
		for(k in 1:N.tag){
				for(i in 1:N.rep){
					lambda[i,1,k,l]		<- betas[1,l] 				+ gammas*COVAR.F[i,1,k,l] + eta[i,1,k,l] + phi[1,k,l] #+ gammas[2]*COVAR.R[i,1,k,l]
				}
  			for(j in 2:N.sp){
				for(i in 1:N.rep){
					lambda[i,j,k,l]		<- betas[1,l] + betas[j,l] 	+ gammas*COVAR.F[i,j,k,l] + eta[i,j,k,l] + phi[j,k,l] #+ gammas[2]*COVAR.R[i,j,k,l] + eta[i,j,k,l] + phi[i,j,k,l]
				}
			}
		}
	}

	for(l in 1:N.t){
		for(k in 1:N.tag){
	  		for(j in 1:N.sp){
	  			for(i in 1:N.rep){
  					eta[i,j,k,l] ~ dnorm(0, 1/sigma2)
				}
 			}
 		}
	}
	for(l in 1:N.t){
 		for(k in 1:N.tag){
	 		for(j in 1:N.sp){
 					phi[j,k,l]	~ dnorm(0, 1/tau2)
			}
		}
	}

	### Derived Quantities
	for(l in 1:N.t){
		p[1,l] <- exp(betas[1,l])
		for(j in 2:N.sp){
			p[j,l]	<-	exp(betas[1,l] + betas[j,l])
		}
	}
	for(l in 1:N.t){
		for(j in 1:N.sp){
			P[j,l]	<-	p[j,l] / sum(p[,l])
		}
  	}

  	### Priors
	for(l in 1:N.t){
	  	for(j in 1:N.sp){
		 	betas[j,l] ~ dnorm(0,0.001)
    	}
	}
#   	for(i in 1:N.gammas){
		 gammas ~ dnorm(0,0.001)
#    	}


 	sigma2 	~ dunif(0,1000)
 	tau2	~ dunif(0,1000)

}
# ",
file="jags_dummy1.txt")

  jags.data     = list("DATA","N.t","N.sp","N.tag","N.rep","COVAR.F") #,"COVAR.R","N.gammas"
  jags.params   = c("betas","gammas","sigma2","tau2","P") #
  model.loc		= c("jags_dummy1.txt")

	Nburn = 100000
	Niter = 10000
	Nchain=	3
	Inits = NULL
	for(i in 1:Nchain){
		Inits[[i]]	<- list("betas"=matrix(runif(N.t*N.sp,-10,10),N.sp,N.t),
 						"gammas"=runif(N.gammas,-10,10),
						"sigma2"=runif(1,0,30),
 						"tau2"=runif(1,0,30))
	}
  	jags.model.rand.v3.b = jags(jags.data, inits = Inits, parameters.to.save= jags.params, model.file=model.loc,
  						n.chains = Nchain, n.burnin = Nburn, n.thin = 1, n.iter = Nburn+Niter, DIC = TRUE)



























# 999999999999999999999999999999999999999999999999999999999999999999999999999999
##########################################################################################
##########################################################################################
## V4
##### FIXED for the match between the tags + primer and the substrate
#####	+ TWO RANDOM EFFECT one shared by all species and one for each time (with realizations for each species and tag)
##########################################################################################

N.gammas	<-	1#length(COVAR)

jagsscript = cat(
# "
model {
	for(l in 1:N.t){
	  	for(k in 1:N.tag){
		  	for(j in 1:N.sp){
	  			for(i in 1:N.rep){
  					DATA[i,j,k,l] ~ dpois(exp(lambda[i,j,k,l]))
				}
			}
 		}
 	}

	for(l in 1:N.t){
		for(k in 1:N.tag){
				for(i in 1:N.rep){
					lambda[i,1,k,l]		<- betas[1,l] 				+ gammas*COVAR.F[i,1,k,l] + eta[i,1,k,l] + phi[1,k,l] #+ gammas[2]*COVAR.R[i,1,k,l]
				}
  			for(j in 2:N.sp){
				for(i in 1:N.rep){
					lambda[i,j,k,l]		<- betas[1,l] + betas[j,l] 	+ gammas*COVAR.F[i,j,k,l] + eta[i,j,k,l] + phi[j,k,l] #+ gammas[2]*COVAR.R[i,j,k,l] + eta[i,j,k,l] + phi[i,j,k,l]
				}
			}
		}
	}

	for(l in 1:N.t){
		for(k in 1:N.tag){
	  		for(j in 1:N.sp){
	  			for(i in 1:N.rep){
  					eta[i,j,k,l] ~ dnorm(0, 1/sigma2)
				}
 			}
 		}
	}


	for(l in 1:N.t){
 		for(k in 1:N.tag){
	 		for(j in 1:N.sp){
 					phi[j,k,l]	~ dnorm(0, 1/tau2[l])
			}
		}
	}



	### Derived Quantities
	for(l in 1:N.t){
		p[1,l] <- exp(betas[1,l])
		for(j in 2:N.sp){
			p[j,l]	<-	exp(betas[1,l] + betas[j,l])
		}
	}
	for(l in 1:N.t){
		for(j in 1:N.sp){
			P[j,l]	<-	p[j,l] / sum(p[,l])
		}
  	}

  	### Priors
	for(l in 1:N.t){
	  	for(j in 1:N.sp){
		 	betas[j,l] ~ dnorm(0,0.001)
    	}
	}
#   	for(i in 1:N.gammas){
		 gammas ~ dnorm(0,0.001)
#    	}


 	sigma2 ~ dunif(0,1000)

  	for(l in 1:N.t){
	 	 	tau2[l] ~ dunif(0,1000)
 	}


}
# ",
file="jags_dummy1.txt")

jags.data     = list("DATA","N.t","N.sp","N.tag","N.rep","COVAR.F") #,"COVAR.R","N.gammas"
jags.params   = c("betas","gammas","sigma2","tau2","P") #
model.loc		= c("jags_dummy1.txt")

Nburn = 100000
Niter = 10000
Nchain=	3
Inits = NULL
for(i in 1:Nchain){
	Inits[[i]]	<- list("betas"=matrix(runif(N.t*N.sp,-10,10),N.sp,N.t),
						"gammas"=runif(N.gammas,-10,10),
						"sigma2"=runif(1,0,30),
						"tau2"=runif(N.t,0,30))
}
jags.model.rand.v4 = jags(jags.data, inits = Inits, parameters.to.save= jags.params, model.file=model.loc, n.chains = Nchain, n.burnin = Nburn, n.thin = 1, n.iter = Nburn+Niter, DIC = TRUE)

































# 10_10_10_10_10_10_10_10_10_10_10_10_10_10_10_10_10_10_10_10_10_10_10_10_10_10_
##########################################################################################
##########################################################################################
## V5
##### FIXED for the match between the tags + primer and the substrate
#####	+ TWO RANDOM EFFECT one shared by all species and one for each time (with realizations for each species and tag)
##########################################################################################

N.gammas	<-	1#length(COVAR)

jagsscript = cat(
# "
model {
	for(l in 1:N.t){
	  	for(k in 1:N.tag){
		  	for(j in 1:N.sp){
	  			for(i in 1:N.rep){
  					DATA[i,j,k,l] ~ dpois(exp(lambda[i,j,k,l]))
				}
			}
 		}
 	}

	for(l in 1:N.t){
		for(k in 1:N.tag){
				for(i in 1:N.rep){
					lambda[i,1,k,l]		<- betas[1,l] 				+ gammas*COVAR.F[i,1,k,l] + eta[i,1,k,l] + nu[1,k,l] #+ phi[1,k,l]  #+ gammas[2]*COVAR.R[i,1,k,l]
				}
  			for(j in 2:N.sp){
				for(i in 1:N.rep){
					lambda[i,j,k,l]		<- betas[1,l] + betas[j,l] 	+ gammas*COVAR.F[i,j,k,l] + eta[i,j,k,l] + nu[j,k,l] #+ phi[j,k,l]  #+ gammas[2]*COVAR.R[i,j,k,l] + eta[i,j,k,l] + phi[i,j,k,l]
				}
			}
		}
	}

	for(l in 1:N.t){
		for(k in 1:N.tag){
	  		for(j in 1:N.sp){
	  			for(i in 1:N.rep){
  					eta[i,j,k,l] ~ dnorm(0, 1/sigma2)
				}
 			}
 		}
	}


# 	for(l in 1:N.t){
#  		for(k in 1:N.tag){
# 	 		for(j in 1:N.sp){
#  					phi[j,k,l]	~ dnorm(0, 1/tau2[l])
# 			}
# 		}
# 	}

	for(l in 1:N.t){
 		for(k in 1:N.tag){
	 		for(j in 1:N.sp){
 					nu[j,k,l]	~ dnorm(0, 1/sigma2.nu[k,l])
			}
		}
	}

	### Derived Quantities
	for(l in 1:N.t){
		p[1,l] <- exp(betas[1,l])
		for(j in 2:N.sp){
			p[j,l]	<-	exp(betas[1,l] + betas[j,l])
		}
	}
	for(l in 1:N.t){
		for(j in 1:N.sp){
			P[j,l]	<-	p[j,l] / sum(p[,l])
		}
  	}

  	### Priors
	for(l in 1:N.t){
	  	for(j in 1:N.sp){
		 	betas[j,l] ~ dnorm(0,0.001)
    	}
	}
#   	for(i in 1:N.gammas){
		 gammas ~ dnorm(0,0.001)
#    	}

	# Variance Priors
 	sigma2 ~ dunif(0,1000)
#   	for(l in 1:N.t){
# 	 	 	tau2[l] ~ dunif(0,1000)
#  	}
	for(l in 1:N.t){
 		for(k in 1:N.tag){
			sigma2.nu[k,l] ~ dunif(0,1000)
		}
	}



}
# ",
file="jags_dummy1.txt")

jags.data     = list("DATA","N.t","N.sp","N.tag","N.rep","COVAR.F") #,"COVAR.R","N.gammas"
jags.params   = c("betas","gammas","sigma2","sigma2.nu","P") #"tau2"
model.loc		= c("jags_dummy1.txt")

Nburn = 100000
Niter = 10000
Nchain=	3
Inits = NULL
for(i in 1:Nchain){
	Inits[[i]]	<- list("betas"=matrix(runif(N.t*N.sp,-10,10),N.sp,N.t),
						"gammas"=runif(N.gammas,-10,10),
						"sigma2"=runif(1,0,30),
#  						"tau2"=runif(N.t,0,30),
						"sigma2.nu"=matrix(runif(N.t*N.tag,0,30),N.tag,N.t))
}
jags.model.rand.v5 = jags(jags.data, inits = Inits, parameters.to.save= jags.params, model.file=model.loc, n.chains = Nchain, n.burnin = Nburn, n.thin = 1, n.iter = Nburn+Niter, DIC = TRUE)






##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

### EXTRACTING INFORMATION FROM THE MODEL OBJECTs













jags.model.rand.base
jags.model.rand.base.2fix
jags.model.rand.base.both
jags.model.base.c

jags.model.rand.v1
jags.model.rand.v2
jags.model.rand.v3
jags.model.rand.v3.b
jags.model.rand.v4
jags.model.rand.v5
