#change wd
setwd(getwd())

#read in the data with allele frequencies
dataChr<-read.table("dataChrPost.txt", header=TRUE, row.names=NULL)

#read in the households data simulated previously to host the initial infections
num.homesteads<-1800
homesteads <- read.table("homesteads.txt", header=TRUE) 

#number of SNPs
numSNPsPerInf<-53

#number of people per household
maxNumPplehh<-8

#load required functions
#get baseline SNPS
getBaselineSNPs<-function(domAlleleFreq) {
  rand<-runif(1, 0.0, 1.0)
  snp<-1
  if (rand>domAlleleFreq) snp<-0
  return(snp)
}

# create a grid to simulate the initial infections
xh <- rep(seq(1:50),50)
yh <- rep(seq(1:50),each=50)
num.loc <- length(xh)

#create grid with columns equal to the number of SNPs 
infections<-array(-9, dim=c(length(xh), numSNPsPerInf))

# get the allele frequancies for the initial infection starting in the bottom LH corner at (1,1) for sweep 1
n.cols <- numSNPsPerInf
domAlleleFreqList<-dataChr$domAlleleFreq
firstInf<-rep(-9,numSNPsPerInf)
for (j in 1:n.cols) {
  firstInf[j]<-getBaselineSNPs(domAlleleFreqList[j])
}

initInf<-infections

#get the allele frequencies for all the other infections in the grid
for (i in 1:nrow(initInf)){
  for (j in 1:n.cols) {
    initInf[i,j]<-getBaselineSNPs(domAlleleFreqList[j])
  }
}

initInf<-cbind(xh,yh,initInf)

# sort the infections by distance 
distI <- sqrt( (xh-50)^2 + (yh-50)^2 )
initInf <- initInf[order(distI),]

#check the scale in the initInf data
library(scales)
initInf[,"xh"]<-rescale(initInf[,"xh"], to = c(0.00,9.00))
initInf[,"yh"]<-rescale(initInf[,"yh"], to = c(0.00,9.00))

# overlay the initial infections with the households
# plot(initInf[,"xh"], initInf[,"yh"],  pch=21, col=as.factor(brks), cex=1.2)
# points(homesteads$locx, homesteads$locy, col="red", pch=16)

# randomly select houses and take 'nearest' infection
numHousesWithInitInf<-13000
infHouses<-sample(seq(1:num.homesteads),numHousesWithInitInf, replace=T)
infHousesInitInf<-array(-9,dim=c(numHousesWithInitInf, numSNPsPerInf+4))
maxDiff<-1.5
for (hh in 1:length(infHouses)) {
    hlocx<-homesteads$locx[infHouses[hh]]
    hlocy<-homesteads$locy[infHouses[hh]]
    hpID<-0
    radiusSet <- initInf[initInf[,"xh"]>(hlocx-maxDiff) & initInf[,"xh"]<(hlocx+maxDiff)
                   & initInf[,"yh"]>(hlocy-maxDiff) & initInf[,"yh"]<(hlocy+maxDiff),]
    if (is.vector(radiusSet)) {
       infH<-radiusSet[3:(numSNPsPerInf+2)]
    }
    if (!is.vector(radiusSet)) {
       sampleI<-sample(1:length(radiusSet[,1]),1)
       infH <-radiusSet[sampleI,3:(numSNPsPerInf+2)]
    }
    infHousesInitInf[hh,1:(numSNPsPerInf+4)]<-c(infHouses[hh],hlocx,hlocy,hpID,infH)
}

infHousesInitInf<-data.frame(infHousesInitInf)
colnames(infHousesInitInf)[1:4]<-c("location.id","locx","locy", "pID")

infHousesInitInf$pID<-sample(1:maxNumPplehh, length(infHousesInitInf[,1]), replace=T)

write.table(infHousesInitInf, "infHousesInitInf.txt",row.names=FALSE)



 
