#species occupancy modeling
#what's the probability of seeing some number of counts of eDNA when the species is not, in fact, present?  i.e., what's the false-positive rate?
#see Lahoz-Monfort 2015 Mol Ecol Res

setwd("/Users/rpk/GoogleDrive/Kelly_Lab/Projects/Carkeek_eDNA_grid/Analysis/Ryan")


### Required packages
# install.packages("unmarked")
library(unmarked)  # Maximum likelihood estimation for SODMs 
# install.packages("jagsUI")
library(jagsUI)    # Bayesian analysis with JAGS
library(data.table)
# This function is later used for printing simulation results on the screen
printsummres<-function(thetahat,thename="estimated parameter"){
  cat(paste0("\n",thename,": mean = ",round(mean(thetahat),2),
              ", 2.5% = ",round(quantile(thetahat,prob=0.025),2),
              ", 97.5% = ",round(quantile(thetahat,prob=0.975),2)))
}

	OTUs=read.csv("/Users/rpk/GoogleDrive/Kelly_Lab/Projects/Carkeek_eDNA_grid/Data/OTU_table.csv", row.names=1)
	meta=read.csv("/Users/rpk/GoogleDrive/Kelly_Lab/Projects/PS_urban_eDNA/AnalysisFinal/metadata_field_exp.csv")
	spatial=meta[meta$experiment=="spatial",]  #select relevant samples from the metadata
	spatial_tags=paste0("lib_", spatial $library,"_tag_", spatial $tag_sequence)  #generate tag-ids for samples from relevant experiment
	all_tags=paste0("lib_", meta$library,"_tag_", meta$tag_sequence)
	ReplicatedExpts=meta[meta$experiment%in%c("control", "spatial"),]
	replicated_tags=paste0("lib_", ReplicatedExpts $library,"_tag_", ReplicatedExpts $tag_sequence)


replicatedOTUs=OTUs[, replicated_tags]
	replicatedOTUs= replicatedOTUs[rowSums(replicatedOTUs)>0,]
	colnames(replicatedOTUs)= ReplicatedExpts$sample_name
	table(droplevels(as.factor(ReplicatedExpts$sample_name)))

##############################################
##############################################
##############################################
##############################################

#needs A binary detection matrix where rows=sites, & columns = replicates.

myRealData = replicatedOTUs ; myRealData[myRealData>0]<-1 ; myRealData=t(myRealData)
#drop NTC because only replicated twice
#myRealData= myRealData[-which(row.names(myRealData)=="NTC."),]


#treating each transect as its own site ; reshape data for model fitting
transectNames=unique(row.names(myRealData))
nTransects=length(unique(transectNames))
nReplicates=4  #replicates per transect/sample/whatever
OTUlist=list()
for (i in 1:ncol(myRealData)){
	OTU=c(myRealData[,i])
		dim(OTU)<-c(nReplicates, nTransects)  #reshape ; note these params will depend upon the input data in question
	OTUlist[[i]]=t(OTU)
}


##############################################
##############################################
##############################################
##############################################
##############################################

####Bayesian version
### Run first this part once to create the file with the JAGS model
sink("RoyleLink_prior.txt")
cat("model {
    # Priors
    psi ~ dunif(0,1)
    p11 ~ dunif(0.1,1)
    p10 ~ dunif(0.03, p10_max)
    
    # Likelihood 
    for (i in 1:S){
    z[i] ~ dbern(psi)
    p[i] <- z[i]*p11 + (1-z[i])*p10
    for (j in 1:K){
    Y[i,j] ~ dbern(p[i])
    }
    }
    } ",fill=TRUE)
sink()


model_Bayesian <- function(datalist,nOTUs=length(datalist), S= nTransects, K= nReplicates, doprint=TRUE,p10_max=0.1,
                                     ni=4000,nt=2,nc=1,nb=400,myparallel=TRUE) {   
  psihat<-p11hat<-p10hat<-rep(nOTUs)
  modelSummaries<-list()
  for(ii in 1: nOTUs){
    if (doprint) cat("\r", ii, "of", nOTUs,"   ")    
    hh<-datalist[[ii]]
    # fit the model    
    jags.inits <-function()(list(psi=runif(1,0.05,0.95),p11=runif(1, 0.1,1),p10=runif(1,0.03,p10_max)))
    jags.data  <-list(Y=hh,S=S,K=K,p10_max=p10_max)
    jags.params<-c("psi","p11","p10")
    model<-jags(data = jags.data, inits = jags.inits, parameters.to.save= jags.params, 
                model.file= "RoyleLink_prior.txt", n.thin= nt, n.chains= nc, 
                n.iter= ni, n.burnin = nb, parallel=myparallel)  #, working.directory= getwd()
    # extract results (medians of the marginal posteriors)
    psihat[ii] <- model$summary["psi","50%"]
    p11hat[ii] <- model$summary["p11","50%"]
    p10hat[ii] <- model$summary["p10","50%"]    
    modelSummaries[[ii]]<-model$summary
  }
  if (doprint){
    printsummres(psihat,thename="estimated psi")
    printsummres(p11hat,thename="estimated p11")
    printsummres(p10hat,thename="estimated p10")
  }
  # pdf(paste0("/Users/rpk/GoogleDrive/Kelly_Lab/Projects/PS_urban_eDNA/Analysis/BayesianModeling/",format(Sys.time(), "%H_%M_%S"),"_p10hat.pdf"))
  # hist(p10hat)
  # dev.off()
  saveRDS(modelSummaries, paste0(format(Sys.time(), "%H_%M_%S"),"_ModelSummaries.rds"))
	BayesResults<-list(psihat=psihat,p11hat=p11hat,p10hat=p10hat,modelSummaries=modelSummaries)
  return(BayesResults)
}

#grabset=sample(1:length(OTUlist), 100)
focalOTUlist=OTUlist  #select dataset from OTUs for analysis

modelSummaries = model_Bayesian(focalOTUlist)  #  run model and store results

##having estimated the parameters, plug into species occupancy model and see what prob of occurrence is given X detections:
#define probability of occupancy, given parameters, according to binomial distribution and equation in Lahoz-Monfort
ProbOcc=function(x, psi, p11, p10, K){
	(psi*(p11^x)*(1-p11)^(K-x)) / ((psi*(p11^x)*(1-p11)^(K-x))+(((1-psi)*(p10^x))*((1-p10)^(K-x))))
	}

ProbOccOut=NA
for (i in 1:length(modelSummaries[[1]])){
	K = 3 #replicates per site 
	psi = modelSummaries[[1]][i] #prob of site occupancy
	p11 = modelSummaries[[2]][i] #prob of detection, given occupancy (i.e., true positive rate)
	p10 = modelSummaries[[3]][i] #prob of detection, given lack of occupancy (i.e., false positive rate)
	nObs= max(rowSums(focalOTUlist[[i]], na.rm=T)[rowSums(focalOTUlist[[i]], na.rm=T)>0]) #sum(OTUlist[[i]], na.rm=T) #treating all occurrences together #min(rowSums(OTUlist[[i]], na.rm=T)[rowSums(OTUlist[[i]], na.rm=T)>0]) #here, using the min number of detections of that OTU across all transects, all sites... essentially asking if there is any transect for which this OTU is likely a false positive.  There are other ways of doing this.
ProbOccOut[i]=ProbOcc(nObs, psi, p11, p10, K)  
}

ProbOccOut #; hist(ProbOccOut)


#focalOTUlist[[which.min(ProbOccOut)]]
ProbOccOut = data.frame(row.names(replicatedOTUs), ProbOccOut) ; names(ProbOccOut)=c("OTUname", "ProbOccupancy")
write.csv(ProbOccOut, paste0(format(Sys.time(), "%H_%M_%S"),"_OTU_probs_out.csv"), row.names=F)


b=ProbOccOut[ProbOccOut[,2]>0.8,1]
d=row.names(OTUs)%in%b
write.csv(OTUs[d,], "OTUs_BayesianVetted_OTUs.csv", row.names=T, quote=F)

