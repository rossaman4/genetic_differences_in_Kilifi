#data preparation
setwd(getwd())

#simulate households and calculate the distance differences
source("simulate_household_distances.R")

#simulate initial infections
source("simulate_initial_infections.R")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# combines functions and simulation to make one R script for sciCORE cluster
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD all the required functions
#_____________________________________________________________________________________________________________
##############  add genoID   #######################
# add genoID identifier (made from individual SNPs)
#_____________________________________________________________________________________________________________


addGenoID<-function(xx) {
  
  genoID<-xx[numSNPsPerInf+numPreSNPcol]
  for (g in 1:(numSNPsPerInf-1)) {
    genoID<-genoID + ( (2^g)*xx[(numSNPsPerInf+numPreSNPcol-g)] )
  }
  xx<-cbind(xx,genoID)
  colnames(xx)[(numSNPsPerInf+numPreSNPcol+1)]<-"genoID"
  
  return(xx)
}

#________________________________________________________________________________________________________________
##############  Function to check  new infections against current and recent ##########################################################
#_______________and delete if necessary_______________________________________________________________________________

checkNewInfections<-function(newInf) {
  
  # newInf<-data.frame(newInf)
  numberNew<-length(newInf[,1])
  if (numberNew>0) {
    # check recentInfections genotypes per house per new infection
    #check recent infs genotypes per house, per person per new infection
    for (nNew in 1:numberNew) {
      # narrow down to the same location
      #and the same person
      recentMatch<-recentInfections[((recentInfections[,"location.id"]==newInf[nNew,"location.id"]) & (recentInfections[,"pID"]==newInf[nNew,"pID"])),]
      # get number of SNPs diff 
      if (length(recentMatch[,1])>0) {
        recentSNPsdiff<-rep(0,length(recentMatch[,1]))
        
        for (r in 1:length(recentSNPsdiff)) {  
          Sdiff<-rep(0,numSNPsPerInf)
          wGenoID1 <- newInf[nNew, "genoID"]
          wGenoID2 <- recentMatch[r,"genoID"]
          for (mn in (numSNPsPerInf):1) {
            infInd1 <- wGenoID1 %% 2
            wGenoID1 <- trunc(wGenoID1/2)       
            infInd2 <- wGenoID2 %% 2
            wGenoID2 <- trunc(wGenoID2/2) 
            if (infInd1 != infInd2) {Sdiff[mn]<-Sdiff[mn]+1}
          }
          recentSNPsdiff[r]<-sum(Sdiff)
        }
        
        # narrow down to SNPS diff greater than specified
        recentMatch<-recentMatch[recentSNPsdiff<=recentSNPsdiffAllowed,]
        if ( length(recentMatch[,1])>numRecentInfAllowed )  {
          newInf[nNew,"clear"]<-1
        }
      } 
    }
    
    # check total number of current infections per house
    for (n3 in 1:length(newInf[,1])) {
      currentMatch<-current_genotypes[((current_genotypes[,"location.id"]==newInf[n3,"location.id"]) & (current_genotypes[,"pID"]==newInf[n3,"pID"])),]
      if ( length(currentMatch[,1])>numCurrentInfAllowed )  {
        newInf[n3,"clear"]<-1
      }
    }
    
  } # closes if numberNew>0
  
  # remove marked new infections 
  if (is.array(newInf)) {
    newInf<-newInf[newInf[,"clear"]==0,]
  }
  
  return(newInf)
}

#_____________________________________________________________________________________________________________
#  Computediff: Function to compute Time and SNP differences 
#_____________________________________________________________________________________________________________
Computediff<-function(Hospdata) {
  for(i in 1:nrow(Hospdata)){
    # cat(i)
    for(j in 1:nrow(Hospdata)){
      if (i<j){
        pair1<-Hospdata[i,"location.id"];pair2<-Hospdata[j,"location.id"]
        geno1<-Hospdata[i,"genoID"]; geno2<-Hospdata[j,"genoID"]
        SNP.diff<-sum(Hospdata[i,firstSNPcol:(firstSNPcol+numSNPsPerInf)]!=Hospdata[j,firstSNPcol:(firstSNPcol+numSNPsPerInf)])
        time.diff<-abs(as.numeric(Hospdata[i,"time"])-as.numeric(Hospdata[j,"time"]))
        start.time<-as.numeric(Hospdata[i,"time"])
        diff1<-cbind(pair1,pair2,SNP.diff,time.diff,start.time,geno1,geno2)
        diffa<-rbind(diffa,diff1)
      }
    }
  }
  return(diffa)
}


#________________________________________________________________________________________________________________
##############  Function to compute Distance between homesteads  #######################################
#________________________________________________________________________________________________________________
distancediff<-function(mydata) {
  diffa<-array(-9,dim=c(1,3))
  for(i in 1:nrow(mydata)){
    #cat(i)
    for(j in 1:nrow(mydata)){
      #if (i<=j){     
        pair1<-i;pair2<-j
        dist.diff<-sqrt((as.numeric(mydata[i,"locy"])-as.numeric(mydata[j,"locy"]))^2 +(as.numeric(mydata[i,"locx"])-as.numeric(mydata[j,"locx"]))^2 )
        diff1<-cbind(pair1,pair2,dist.diff)  #
        diffa<-rbind(diffa,diff1)
      #}
    }
  }
  diffa<-diffa[diffa[,"pair1"]!=-9,]
  return(diffa)
}

# ___________________________________________________________________
#
# Baseline SNP

getBaselineSNPs<-function(domAlleleFreq) {
  rand<-runif(1)
  snp<-1
  if (rand>domAlleleFreq) snp<-0
  return(snp)
}


#________________________________________________________________________________________________________________
##############  Function to Sample hospital cases  #############################################################
#________________________________________________________________________________________________________________

# get infections to match observeds
getHospSamples<-function(current_genotypes){
  # get distances to house with observed infection
  locx<-obsInf$locx[obst];  locy<-obsInf$locy[obst]; set<-obsInf$set[obst]; setInf<-obsInf$setInf[obst]
  obsInfloc<-data.frame(locx,locy, set, setInf); currGenoloc<-data.frame(current_genotypes$locx, current_genotypes$locy)
  colnames(currGenoloc)<-c("locx", "locy")
  
  dist_genotypes<-current_genotypes
  dist_genotypes$dist<-sqrt((as.numeric(obsInfloc$locy)-as.numeric(currGenoloc$locy))^2 +(as.numeric(obsInfloc$locx)-as.numeric(currGenoloc$locx))^2 )
  dist_genotypes$set<-obsInfloc$set; dist_genotypes$setInf<-obsInfloc$setInf
  
 
  # select closest 100 infections for easy management and then account for infWeight
  dist_genotypes<-dist_genotypes[order(dist_genotypes$dist),]
  dist_genotypes<-dist_genotypes[1:100,]

  # avoid oversampling sibling infections  
  for (dg in 1:length(dist_genotypes[,1])) {
    if (dist_genotypes$infWeight[dg]==0.25) {
      temp3<-dist_genotypes[dist_genotypes$pID==dist_genotypes$pID[dg] & dist_genotypes$location.id==dist_genotypes$location.id[dg] &
                              dist_genotypes$infWeight==0.25 & dist_genotypes$t==dist_genotypes$t[dg],]
      dist_genotypes<-dist_genotypes[!(dist_genotypes$pID==dist_genotypes$pID[dg] & dist_genotypes$location.id==dist_genotypes$location.id[dg] &
                                       dist_genotypes$infWeight==0.25 & dist_genotypes$t==dist_genotypes$t[dg]),]
      temp3<-sample(temp3,1)
      dist_genotypes<-rbind(temp3, dist_genotypes)      
    }
  }
  # adjustment to avoid oversampling recombination siblings
  #dist_genotypes$rand<-runif(length(dist_genotypes[,1]),0,1)
  #dist_genotypes<-dist_genotypes[((dist_genotypes$infWeight==1) | (dist_genotypes$infWeight==0.25 & dist_genotypes$rand<0.25)),]
  #dist_genotypes$rand<-NULL
  
  # select 10 closest infections (max one per person & 10 simulated per observed)  
  dist_genotypes<-dist_genotypes[!duplicated(dist_genotypes[1:2]),]
  dist_genotypes<-dist_genotypes[order(dist_genotypes$dist),]
  hospsample <- dist_genotypes[1:10,]
  
  
  # label closest infection with obst
  hospsample$numCommInf<-obst
  hospsample$time<-t
  # store in Hospdata
  Hospdata<-rbind(Hospdata,hospsample)
  return(Hospdata)
}


#________________________________________________________________________________________________________________
##############  Function to effect mutation  #############################################################
#________________________________________________________________________________________________________________
getMutation<-function(currentSNP,probMutatePerTransmPerSNP) {
  rand1<-runif(1)
  newSNP<-currentSNP   
  if (rand1<probMutatePerTransmPerSNP) {
    possSNPvalues<-c(0,1)
    possSNPvalues<-possSNPvalues[possSNPvalues!=currentSNP]
    newSNP<-possSNPvalues[1]
  } 
  return(newSNP)
}

#________________________________________________________________________________________________________________
##############  Function to get new infections per timestep  ###############################################################
#________________________________________________________________________________________________________________
getNewInfections<-function(t) {
  
   currentNumInf<-length(current_genotypes[,1])

  # dummy row for newInf, omit later after loop
  newInf<-current_genotypes[1,]  
  
  for (i in 1:currentNumInf){
    
    # number of new infections for current infection
    meanNew<-meanNumNewInfPerTimestep*current_genotypes[i,"infWeight"]

    #5% for the infections at the edge of the grid
    if ((current_genotypes[i, "locx"]<0.45) | (current_genotypes[i, "locx"]>8.55) | (current_genotypes[i, "locy"]<0.45) | (current_genotypes[i, "locy"]>8.55)){
      meanNew<-(meanNumNewInfPerTimestep/2)
    }
    
    currenttime<-current_genotypes[i,"time"]
    ageInf<- t - currenttime

    currenthome<-current_genotypes[i,"location.id"]
    currentperson<-current_genotypes[i,"pID"]

    #infectivity over time - decays with an exponential function (an individual is fully infectious for "infsTime" days)
    infsTime<-75
    rho<-1/infsTime
    
    meanNew[currenthome]<-meanNew*exp(-rho*ageInf)
    meanNew[currenthome]<-meanNew[currenthome]*(dist_kernel[currenthome,1])/max(dist_kernel$sum_kernel)
      
    numNewInfPerTimestep<-rpois(1,meanNew[currenthome])
      
    # assign new infections
    if (numNewInfPerTimestep>0){
      for (j in 1:numNewInfPerTimestep) {
        newInf<-rbind(newInf,current_genotypes[i,])
        newInf[length(newInf[,1]),"infWeight"]<-1
        
        # get new locations
     
        newInf[length(newInf[,1]), "location.id"]<-getNewLoc(currenthome)
        
        #condition for person ID - if it is in the same house at the same time then a different person is infected

        maxNumpplehh<-8

	      numPple<-1:maxNumpplehh
        
        if (newInf[length(newInf[,1]), "location.id"]==currenthome & newInf[length(newInf[,1]), "time"]==currenttime){
                numPple<-subset(numPple, numPple!=currentperson)}

		newInf[length(newInf[,1]), "locx"]<-homesteads$locx[homesteads$id==newInf[length(newInf[,1]), "location.id"]]
		newInf[length(newInf[,1]), "locy"]<-homesteads$locy[homesteads$id==newInf[length(newInf[,1]), "location.id"]]

	if (newInf[length(newInf[,1]), "location.id"]!=currenthome & newInf[length(newInf[,1]), "time"]==currenttime){
                numPple<-1:maxNumpplehh}

		newInf[length(newInf[,1]), "locx"]<-homesteads$locx[homesteads$id==newInf[length(newInf[,1]), "location.id"]]
		newInf[length(newInf[,1]), "locy"]<-homesteads$locy[homesteads$id==newInf[length(newInf[,1]), "location.id"]]
	    	
	  newInf[length(newInf[,1]), "pID"]<-sample(numPple,1)


        # get mutations
        for (k in 1:numSNPsPerInf){
          newInf[length(newInf[,1]),numPreSNPcol+k]<-getMutation(newInf[length(newInf[,1]),numPreSNPcol+k],probMutatePerSNPPert)
        }              
        
        #add recombination
        rand3<-runif(1)
        if (rand3<probRecombination){  
          potRecombInf<-current_genotypes[which(current_genotypes[,"location.id"]==currenthome & current_genotypes[,"pID"]==currentperson),]
          if (is.data.frame(potRecombInf)) {
            idR<-seq(1:length(potRecombInf[,1])) # random genotype to recombine with
            idR <- sample(idR, 1,replace=FALSE, prob=potRecombInf[,"infWeight"])  
            rndmgenotype<-potRecombInf[idR,]
           
            # store a further 3 genotypes to newInf to give space for the four siblings and change infWeight to 0.25
            tpRow<-length(newInf[,1])             
            newInf[tpRow,"infWeight"]<-0.25
            newInf<-rbind(newInf, newInf[tpRow,], newInf[tpRow,], newInf[tpRow,])         
            tpRow<-tpRow-1
                 
             for (h in 1:numChrPerInf){ 
               # order 1,2,3,4 at random for each chromosome
               inf_order<-sample(c(1,2,3,4))

               # chromosome from recombining infection from potRecombInf
               newInf[(tpRow+inf_order[1]),c((numPreSNPcol+startChrcol[h]):(numPreSNPcol+endChrcol[h]))]<- rndmgenotype[c((numPreSNPcol+startChrcol[h]):(numPreSNPcol+endChrcol[h]))]        
               
               # chromosome - recombined		 
               newInf[(tpRow+inf_order[2]),c((numPreSNPcol+startChrcol[h]):(numPreSNPcol+endChrcol[h]))]<-getRecombinant(newInf[(tpRow+inf_order[2]), 
		   c((numPreSNPcol+startChrcol[h]):(numPreSNPcol+endChrcol[h]))],rndmgenotype[c((numPreSNPcol+startChrcol[h]):(numPreSNPcol+endChrcol[h]))], 
		   dataChr$cumProbG[c(startChrcol[h]:endChrcol[h])])      
                        
               # third is current chromosome (no change)
               
               # fourth is converse of recombined genotype
              for (fc in 1:(endChrcol[h]-startChrcol[h]+1)) {
                 if (newInf[(tpRow+inf_order[2]),(numPreSNPcol+startChrcol[h]+fc-1)]==newInf[(tpRow+inf_order[3]),(numPreSNPcol+startChrcol[h]+fc-1)] ){
                         newInf[(tpRow+inf_order[4]),(numPreSNPcol+startChrcol[h]+fc-1)]<-rndmgenotype[(numPreSNPcol+startChrcol[h]+fc-1)]
                 }         
              }
             }                           
           }
        }

        # get new genoID (code based on addGenoID)
        genoID<-newInf[,(numSNPsPerInf+numPreSNPcol)]
        for (g in 1:(numSNPsPerInf-1)) {
          genoID<-genoID + ( (2^g)*newInf[,(numSNPsPerInf+numPreSNPcol-g)] )
        }
        newInf[,"genoID"]<-genoID
        
        # set clearance to 0 for new infections
        newInf[,"clear"]<-0 
        # get timestep
        newInf[,"time"]<-rep(t,length(newInf[,1]))
        # col1=clear, col2=time, col3=locx, col4=locy, col5-104=snps 
      }
    }
  } # closes i loop
  
  # tidy new infections up
  # remove first dummy row of newInf
  if (is.array(newInf)) {
    newInf<-newInf[-1,]
  }  
  
  return(newInf)
} 


#________________________________________________________________________________________________________________
##############  Function to  select new household for the new infection ###############################################################
#________________________________________________________________________________________________________________
getNewLoc<-function(currenthome) {
  specificHomes <- dist_diff[ which(dist_diff$id==currenthome | dist_diff$id_x ==currenthome), ]
  specificHomes[,"scaledProp"]<-specificHomes[,"prop"]/sum(specificHomes[,"prop"])
  specificHomes <- specificHomes[order(-specificHomes[,"scaledProp"]),]
  specificHomes[,"cumProp"]<-cumsum(specificHomes[,"scaledProp"])
  rand_new<-runif(1)
  
  specificHomes<-specificHomes[c(1,1:nrow(specificHomes)),]
  specificHomes[1,"cumProp"]<-0
  
  smallHomes<-specificHomes[specificHomes[,"cumProp"]<rand_new,]
  homePair<-smallHomes[nrow(smallHomes),1:2]
  
  newloc<-homePair[1]
  if (newloc==currenthome) newloc<-homePair[2]
  
  return(newloc)
}


# __________________________________________________________________________
#
# getPotentialImportedInfections
# function to create a set of imported infections which can potentially be added in to the simulated infections

getPotentialImportedInfections<-function() {

  nImp<-2
  
  potAddInf <- array(-9, dim=c(nImp, numSNPsPerInf))
  for (i in 1:nImp) {
    for (j in 1:numSNPsPerInf) {
      potAddInf[i,j]<-getBaselineSNPs(domAlleleFreqList[j]) 
    }
  }
  
  potAdd<-sample(homesteads[,"id"], nImp, replace=TRUE)
  potAdd<-homesteads[potAdd,]
  potAdd$pID<-sample(1:maxNumpplehh, nImp, replace=T)
  colnames(potAdd)<-c("potAddH", "locx", "locy","pID")
  
  # add a list of random times (0-2500) 
  potAddt <- sample(1:2850, nImp, replace=TRUE) #time changed to reflect new timesteps in obsInf
  potAddt <- sort(potAddt)
  # make last time very long to avoid missing value problem
  potAddt[length(potAddt)]<-999999
  clear<- rep(0,nImp)
  numCommInf <- rep(0,nImp)
  infWeight<- rep(1, nImp)
  potAddInf <- cbind(potAdd$potAddH, potAdd$pID, clear, potAddt, numCommInf, potAdd$locx, potAdd$locy, infWeight, potAddInf)
  potAddInf<-data.frame(potAddInf)
  names(potAddInf)[1]<-"potAddH";names(potAddInf)[2]<-"pID";names(potAddInf)[6]<-"locx";names(potAddInf)[7]<-"locy"; names(potAddInf)[8]<-"infWeight"

  return(potAddInf)
}


#____________________________________________________________________________________________________________
##############  Function to effect Recombination  ###############################################################
#____________________________________________________________________________________________________________
getRecombinant<-function(currentChr,randomChr, cumProbG) {
  rand4<-runif(1)
  breakpoint<-0
  
  for(p in 1:length(currentChr)){
    if(rand4<cumProbG[p]){breakpoint<-p}
  }
  
  recoChr<-currentChr
    if(breakpoint>0){
    recoChr[1:breakpoint]<-randomChr[1:breakpoint]
  }
  
  return(recoChr)
}
# _______________________________________________________________________
# setUpInitialCurrentGenotypes
# function to set up initial current_genotypes from read-in initial infections
# ----------------------------------
numPreSNPcol<-8
setUpInitialCurrentGenotypes<-function() {
  num.init<-dim(infHousesInitInf)[1]
  clear<-rep(0,num.init)
  time<-rep(0,num.init)
  numCommInf<-rep(num.init, num.init)
  current_genotypes <- array(-9, dim=c(length(infHousesInitInf[,1]), numSNPsPerInf+numPreSNPcol))
  current_genotypes[,1] <- infHousesInitInf[,"location.id"]
  current_genotypes[,2] <- infHousesInitInf[,"pID"]
  current_genotypes[,3] <- clear
  current_genotypes[,4] <- time
  current_genotypes[,5] <- numCommInf
  current_genotypes[,6] <- infHousesInitInf[,"locx"]
  current_genotypes[,7] <- infHousesInitInf[,"locy"]
  current_genotypes[,8] <- 1
  for (i in 1:numSNPsPerInf) {current_genotypes[,(i+numPreSNPcol)] <- infHousesInitInf[,(i+4)]}
  current_genotypes<-data.frame(current_genotypes)
  colnames(current_genotypes)[1:numPreSNPcol]<-c("location.id","pID","clear","time","numCommInf", "locx", "locy", "infWeight")
  
  firstSNPcol<-numPreSNPcol+1
  current_genotypes<-addGenoID(current_genotypes)
}


#________________________________________________________________________________________________________________
##############   updateCurrentGenotypes for one time step     ##############################################
#________________________________________________________________________________________________________________

updateCurrentGenotypes<-function(t){
  
  # mark for clearance
  for (md in 1:length(current_genotypes[,1])) {    
    rand<-runif(1)
    if (rand<probClearancePerTimestep){current_genotypes[md,"clear"]<-1}
    if (t-current_genotypes[md,"time"]>maxDuration) {current_genotypes[md,"clear"]<-1}    
  }
  # for recombined infections with siblings, need to clear as a group 
  for (md in 1:length(current_genotypes[,1])) {  
     if ((current_genotypes[md,"clear"]==1) & (current_genotypes[md,"infWeight"]<1)) {
        temppID<-current_genotypes[md,"pID"]
        tempInfWeight<-current_genotypes[md,"infWeight"]
        tempTime<-current_genotypes[md,"time"]
        tempLoc<-current_genotypes[md,"location.id"]
        tempSample<-current_genotypes[current_genotypes["pID"]==temppID & current_genotypes["time"]==tempTime & 
                          current_genotypes["location.id"]==tempLoc & current_genotypes["infWeight"]==tempInfWeight,"clear"]
        current_genotypes[current_genotypes["pID"]==temppID & current_genotypes["time"]==tempTime & 
                          current_genotypes["location.id"]==tempLoc & current_genotypes["infWeight"]==tempInfWeight,"clear"]<-sample(tempSample,1,replace=FALSE)                
      }  
  }

  # clear doomed infections 
  current_genotypes<-current_genotypes[current_genotypes[,"clear"]==0,]
  
  # make new "current" dataset
  current_genotypes<-rbind(current_genotypes,newInf)
  
  # to speedup the iteration set ceiling
 # if (length(current_genotypes[,1])>25000){
 #   nt<-seq(1:length(current_genotypes[,1]))
 #   nt <- sample(nt, 20000,replace=FALSE)
 #   current_genotypes<-current_genotypes[nt,]
    #cat("ceiling reached at time",t,"\n")
#  }
  
  return(current_genotypes)
}


#_____________________________________________________________________________________________________________
##############  Function to update recentInfections  
#____________________________________________________________________________________________________________


updateRecentInfections<-function(t) {
  
  # add new infections to recentInfections store
  recentInfections<-rbind(recentInfections,newInf[c(1,2,4,6,7,61)])
  # remove older past infections
  recentInfections<-recentInfections[(recentInfections[,"time"]>(t-recentInfTime)),]
  return(recentInfections)
}
# End of functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Load required libraries
library(SpatialEpi);library(data.table); require(utils)
library(plyr);library(dplyr); library(dtplyr)
library(foreign);library("scales")

############################################################################################\
#Create subfolders to save files

subDir<-"base Model"

dir.create(paste(getwd(),"/", subDir, sep=""))

############################################################################################
#SET UP THE PARAMETER VALUES

rate<-0.5001
# 'rate' is the input parameter for the probability of the location of an offspring infection by distance.
# It the value of sigma in a half-normal distribution
# The mean of the half-normal distribution is given by sigma*(sqrt(2/pi))
# The original publication erroneously described sigma as the mean, but this has been corrected
probRecombination<-0.5001
meanNumNewInfPerTimestepI<-0.0506
meanNumNewInfPerTimestepW<-0.0485
meanNumNewInfPerTimestepR1<-0.0462
meanNumNewInfPerTimestepR2<-0.0414

# OBSERVED INFECTIONS
setwd(getwd())

#read the observed data
obsInf<-read.csv("obsData.csv", header=TRUE, row.names=NULL)

# HOUSE-PAIRS AND DISTANCES 
# allpairs is a dataset with simulated house locations and obsInf locations and distances between every possible pair of house-obsInf
# to find closest house with an infection at time t 
allpairs<-read.table("allpairs_dist.txt", header=TRUE, row.names=NULL)

# HOMESTEADS
# simulate homesteads and calculate distances for infection 'movement', if not previously done and no homestead data
num.homesteads<-1800

# maximum number of people in a household (average)
maxNumpplehh<-8

#read in simulated homesteads
homesteads<-read.table("homesteads.txt", header=T, row.names = NULL)

#read in file for chromosome and SNP positions
dataChr<-read.table("dataChrPost.txt", header=T, row.names = NULL)

numChrPerInf<-length(unique(dataChr$chromosome))
chrlist<-unique(dataChr$chromosome)
startChrcol<-rep(0, numChrPerInf)
endChrcol<-rep(0, numChrPerInf)
for(zz in 1:numChrPerInf){
	currentChr<-chrlist[zz]
	for(w in 1:length(dataChr$chromosome)){
		if(dataChr$chromosome[w]==currentChr){
		endChrcol[zz]<-w
		}
	}
	for(g in length(dataChr$chromosome):1){
		if(dataChr$chromosome[g]==currentChr){
		startChrcol[zz]<-g
		}
	}
}
	
dist_diff<-read.table("dist_diff.txt", header=TRUE, row.names=NULL)


# 'rate' is an input parameter for distance (set above or read in)
mu<-1/(2*(rate^2))
dist_diff[,"prop"]<-exp(-mu*(dist_diff[,"dist.diff"]^2)) #the normal kernel

distlist = list()
for (i in 1:num.homesteads){
  homestes<-dist_diff[ which(dist_diff$id==i | dist_diff$id_x ==i), ]
  distlist[[i]] <- sum(homestes$prop)
}

dist_kernel = do.call(rbind, distlist)
dist_kernel<-data.frame(dist_kernel, seq(1:1800))
names(dist_kernel)<-c("sum_kernel", "house")


# number of SNPs
numSNPsPerInf<-53

# max duration of infections (otherwise exponential) 
maxDuration<-4000

# load the initial infections

infHousesInitInf<-read.table("infHousesInitInf.txt",header=TRUE)
colnames(infHousesInitInf)[1:4]<-c("location.id","locx","locy", "pID")

current_genotypes<-setUpInitialCurrentGenotypes() 
genoID<-current_genotypes[numSNPsPerInf+numPreSNPcol+1]

# store for recent infections (option to avoid re-infecting same house many times with same genotype)
# label is genoID
recentInfections<-cbind(current_genotypes[,c(1,2,4,6,7,8)],genoID)
colnames(recentInfections)<-c("location.id","pID", "time", "locx", "locy", "infWeight","genoID")

# allele frequencies for imported infections
domAlleleFreqList<-dataChr$domAlleleFreq

# set up array of potential imported infections to be added in to the simulation
potAddInf<-getPotentialImportedInfections()
potAddInf<-addGenoID(potAddInf)

# UPDATE TIMESTEPS AND GET HOSPITAL SAMPLES

# store for new infections (delete first row afterwards)
Hospdata<-cbind(current_genotypes[1,], 0, 0, 0)
nColHD<-dim(Hospdata)[2]
names(Hospdata)[63]<-"dist"
names(Hospdata)[64]<-"set"
names(Hospdata)[65]<-"setInf"


# infection characteristics
currentNumInf<-length(current_genotypes[,1])
meanNumNewInfPerTimestep<-meanNumNewInfPerTimestepW
probClearancePerTimestep<-0.025
probMutatePerSNPPert<-0.0000002

obsInf$time<-obsInf$time+1000
pt<-1

# warm-up period to get a substantial number of infections

for (t in 1:999) {
 newInf<-getNewInfections(t) 
 current_genotypes<-updateCurrentGenotypes(t)
if (t==potAddInf$potAddt[pt]) {
    potImport<-potAddInf[pt,]; names(potImport)<-names(current_genotypes)
    current_genotypes<-rbind(current_genotypes,potImport)
    # currently  included to test 
    pt<-pt+1
  }

if (t==999) {
    # write to text file for investigations
    infs<-current_genotypes
  }
}

fName1 <- paste0(getwd(),"/",subDir,"/","stochInf1","_rate", rate, "_rec",probRecombination,".txt")
write.table(infs,file=fName1, row.names=FALSE)

# stable number of infections
meanNumNewInfPerTimestep<-meanNumNewInfPerTimestepR1
probClearancePerTimestep<-0.025

obst<-1
for (t in 1000:1400) {
  newInf<-getNewInfections(t)				
  current_genotypes<-updateCurrentGenotypes(t)
  if (t==potAddInf$potAddt[pt]) {
    potImport<-potAddInf[pt,]; names(potImport)<-names(current_genotypes)
    current_genotypes<-rbind(current_genotypes,potImport)
    # currently  included to test 
    pt<-pt+1
  }
  
  # sample simulated infections to match time and location of observed infections
  if (t==obsInf$time[obst]) {        
    # get number of observed infections at time t 
    numObst<-nrow(obsInf[obsInf$time==t,])
    # for each observed infection, sample a simulated infection
    for (z in 1:numObst) {
      Hospdata<-getHospSamples(current_genotypes)
      obst<-obst+1 
    }    
  }
  
  if (t==1200) {
    # write to text file for investigations
    infs00<-current_genotypes
  }
  if (t==1400) {
    # write to text file for investigations
    infs01<-current_genotypes
  }
}


#############################
#output files

fName <- paste0(getwd(),"/",subDir,"/","simHosp","_rate", rate, "_mut",probMutatePerSNPPert,"_rec",probRecombination,".txt")

fName00<- paste0(getwd(),"/",subDir,"/","simInf00","_rate", rate, "_rec",probRecombination,".txt")
fName01<- paste0(getwd(),"/",subDir,"/","simInf01","_rate", rate, "_rec",probRecombination,".txt")

write.table(Hospdata,file=fName, row.names=FALSE)
write.table(infs00,file=fName00, row.names=FALSE)
write.table(infs01,file=fName01, row.names=FALSE)


# stable number of infections
meanNumNewInfPerTimestep<-meanNumNewInfPerTimestepR2
probClearancePerTimestep<-0.025

for (t in 1401:1850) {
  newInf<-getNewInfections(t)				
  current_genotypes<-updateCurrentGenotypes(t)
  if (t==potAddInf$potAddt[pt]) {
    potImport<-potAddInf[pt,]; names(potImport)<-names(current_genotypes)
    current_genotypes<-rbind(current_genotypes,potImport)
    # currently  included to test 
    pt<-pt+1
  }
  
  # sample simulated infections to match time and location of observed infections
  if (t==obsInf$time[obst]) {        
    # get number of observed infections at time t 
    numObst<-nrow(obsInf[obsInf$time==t,])
    # for each observed infection, sample a simulated infection
    for (z in 1:numObst) {
      Hospdata<-getHospSamples(current_genotypes)
      obst<-obst+1 
    }    
  }
  
  if (t==1650) {
    # write to text file for investigations
    infs02<-current_genotypes
  }
  if (t==1850) {
    # write to text file for investigations
    infs03<-current_genotypes
  }
}



fName02<- paste0(getwd(),"/",subDir,"/","simInf02","_rate", rate, "_rec",probRecombination,".txt")
fName03<- paste0(getwd(),"/",subDir,"/","simInf03","_rate", rate, "_rec",probRecombination,".txt")

fName <- paste0(getwd(),"/",subDir,"/","simHosp","_rate", rate, "_mut",probMutatePerSNPPert,"_rec",probRecombination,".txt")

write.table(Hospdata,file=fName, row.names=FALSE)
write.table(infs02,file=fName02, row.names=FALSE)
write.table(infs03,file=fName03, row.names=FALSE)

