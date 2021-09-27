############################################################################################################
rm(list=ls())
setwd(getwd())

library(data.table);library(plyr);library(dplyr);library(foreign)
library("scales");library(sampling)

###
numSNPsPerInf=53
getPairs<-function(simData) {
  sim.pairs<-rep(-9,10)
  numSets<-max(simData$set)
  for (i in 1:numSets) {
    simDataSet<-simData[simData$set==i,]
    numInSet<-dim(simDataSet)[1]
    pairs1<-do.call(rbind, replicate(numInSet,simDataSet,simplify=FALSE))
    pairs2<-data.frame(pairs1[order(pairs1$setInf),])
    colnames(pairs1)<-c("set1","setInf1","genoID1")
    colnames(pairs2)<-c("set2","setInf2","genoID2")
    pairs<-cbind(pairs1, pairs2)
    pairs<-pairs[pairs$setInf1 < pairs$setInf2,]
    sim.pairs<-rbind(sim.pairs,pairs)
  }
  sim.pairs<-sim.pairs[sim.pairs$set1!=-9,]
  
  # get SNP differences for data
  SNPsdiff<-rep(0,length(sim.pairs[,1]))
  
  for (i in 1:length(SNPsdiff)) {  
    Sdiff<-rep(0,numSNPsPerInf)
    wGenoID1 <- sim.pairs$genoID1[i]
    wGenoID2 <- sim.pairs$genoID2[i]
    for (mn in (numSNPsPerInf):1) {
      infInd1 <- wGenoID1 %% 2
      wGenoID1 <- trunc(wGenoID1/2)       
      infInd2 <- wGenoID2 %% 2
      wGenoID2 <- trunc(wGenoID2/2) 
      if (infInd1 != infInd2) {Sdiff[mn]<-Sdiff[mn]+1}
    }
    SNPsdiff[i]<-sum(Sdiff)
  }
  sim.pairs<-cbind(sim.pairs, SNPsdiff)
  return(sim.pairs)
}

########################################################################################################
# read in the observed csv files and pair all observations by set
obs.data<-read.csv("obsData.csv", header=TRUE, row.names=NULL)

# get paired observations
obsdata<-subset(obs.data, select = c(set,setInf,genoID))
obspairs<-getPairs(obsdata)

obs.pairs1<-obspairs[,c("set1", "setInf1", "set2","setInf2","SNPsdiff")]

#########################################################################################################
#create a grid for the loglikelihood
para_grid<-data.frame(matrix(0, nrow=1, ncol=4))
para_grid[,1]<-0.5001
para_grid[,2]<-0.5001
names(para_grid)[1:4]<-c("rateList", "recombList","logsum","ssesum")

# read the simulated files by either specific value for distance
# THIS IS A SPECIFIC EXAMPLE ('rate' value of 0.5001)
fileS<-list.files(path= paste0(getwd(), "/",sep=""), pattern='*rate0.5001*')
fileS<-fileS[order(nchar(fileS), fileS)]

obs.pairs<-arrange(obs.pairs1, set1, setInf1, set2, setInf2)
names(obs.pairs)[5]<-"SNPsdiff_obs"

rep<-2
for (i in 1:length(fileS)){
  sim.hosp<-lapply(fileS[i], function(x)read.table(x, header=T))
  simHosp <- ldply(sim.hosp, data.frame)
  simHosp<-simHosp[-1,]
  setDT(simHosp)[, repInfs := seq_len(.N), by = numCommInf]
  simHosp<-data.frame(subset(simHosp,simHosp$repInfs<=rep))

  for (k in 1:rep){
    simhospi<-data.frame(subset(simHosp,simHosp$repInfs==k))
    sim.data<-subset(simhospi, select = c(set,setInf,genoID))
    sim.pairs<-getPairs(sim.data)

    sim.pairs<-sim.pairs[,c("set1", "setInf1", "set2","setInf2","SNPsdiff")]
    
    obs.pairs<-merge(obs.pairs,sim.pairs, by=c("set1", "setInf1", "set2","setInf2"), all=T)
    obs.pairs<-obs.pairs[complete.cases(obs.pairs),]
    obs.pairs<-arrange(obs.pairs, set1, setInf1, set2,setInf2)
    
  }
  
}

obs.pairs$meanSNP=(rowSums(obs.pairs[,c(6:(ncol(obs.pairs)))]))/(ncol(obs.pairs)-5)
obs.pairs$meanSNP[obs.pairs$meanSNP==0]<-0.01
obs.pairs$SNPsdiff_obs[obs.pairs$SNPsdiff_obs==0]<-0.01
obs.pairs$loglik<-(obs.pairs$SNPsdiff_obs*log(obs.pairs$meanSNP/numSNPsPerInf))+((numSNPsPerInf-obs.pairs$SNPsdiff_obs)* log(1-(obs.pairs$meanSNP/numSNPsPerInf)))
#obs.pairs$loglik<-(obs.pairs$SNPsdiff_obs*log(obs.pairs$meanSNP)) - obs.pairs$meanSNP
obs.pairs$resid<-obs.pairs$SNPsdiff_obs - obs.pairs$meanSNP
obs.pairs$sse<-(obs.pairs$SNPsdiff_obs - obs.pairs$meanSNP)^2

logsum<-sum(obs.pairs$loglik)
ssesum<-sum(obs.pairs$sse)

fname<-fileS[i]
y<-as.numeric(unlist(regmatches(fname,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",fname))))
para_grid[y[1]==para_grid$rateList & y[4]==para_grid$recombList,3:4]<-c(logsum,ssesum)

para_grid<-data.frame(para_grid)

#rate<-0.50

fName<- paste0(getwd(),"/","SNPdiffs","_rate", rate, ".txt")
write.table(obs.pairs,file=fName, row.names=FALSE)

write.csv(para_grid, "para_loglik.txt", row.names = F)

