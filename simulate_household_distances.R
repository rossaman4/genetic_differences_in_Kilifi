setwd(getwd())
# Homesteads in simulation evaluation
# simulate 1800 homesteads on a 9by9 grid

num.homesteads<-1800
homesteads<-array(-9, dim=c(num.homesteads,3))
colnames(homesteads)<-c("id","locx","locy")
homesteads[,1]<-seq(1:num.homesteads)
homesteads[,2]<-runif(num.homesteads,0.00,9.00) 
homesteads[,3]<-runif(num.homesteads,0.00,9.00) 
homesteads<-data.frame(homesteads)

maxNumpplehh<-8
#homesteads$pID<-sample(1:maxNumpplehh, num.homesteads, replace=T)
write.table(homesteads,"homesteads.txt", row.names = F)

# Compute distances between pairs of homesteads
# and read this in as a dataset

#function to calculate distance between pairs in kms
distance_km<-function(i,j,data_km){
  dist<-sqrt((data_km$x[i]-data_km$x[j])^2+(data_km$y[i]-data_km$y[j])^2)
  return(dist)
}

dist_km<-Vectorize(distance_km,vectorize.args=list("i","j"))

#calculate the distance between the simulated homesteads
x<-homesteads$locx
y<-homesteads$locy
data_km<-data.frame(x=x,y=y)
time_km<-system.time(res_km<-outer(1:num.homesteads,1:num.homesteads,dist_km,data_km=data_km))


id<-1:nrow(res_km); res_km<-cbind(res_km, id) 

#to get the distance difference for all the pairs
library(reshape2)
distDiff <- melt(as.data.frame(res_km), id="id")

distDiff$id_x<-as.numeric(distDiff$variable); distDiff$variable<-NULL

houses<-homesteads
distanceDiff<-merge(houses, distDiff, by = "id")

colnames(houses)<-c("id_x", "locx_x", "locy_x")
distanceDiff<-merge(houses, distanceDiff, by = "id_x"); names(distanceDiff)[7]<-"dist.diff"
distanceDiff <- subset(distanceDiff, select=c(id_x,id,dist.diff,locx,locy,locx_x,locy_x))

write.table(distanceDiff,"allpairs_dist.txt", row.names = F)

#to get the distance difference for i<j
#remove half the triangle
library(gdata)
res_km1<-res_km
res_km1[upper.tri(res_km1,diag=FALSE)] <- NA
#reshaping the data into long format
res_km1<-as.data.frame(res_km1)
res_km1$id<-1:nrow(res_km1)

library(reshape)
diffDist<-melt(res_km1,id=c("id"))

dist_diff<-diffDist[complete.cases(diffDist), ]
dist_diff$id_x<-as.numeric(dist_diff$variable); dist_diff$variable<-NULL

houses<-homesteads
dist.Diff<-merge(houses, dist_diff, by = "id")

colnames(houses)<-c("id_x", "locx_x", "locy_x")
dist.Diff<-merge(houses, dist.Diff, by = "id_x"); names(dist.Diff)[7]<-"dist.diff"
dist.Diff <- subset(dist.Diff, select=c(id_x,id,dist.diff,locx,locy,locx_x,locy_x))

write.table(dist.Diff,"dist_diff.txt", row.names = F)



