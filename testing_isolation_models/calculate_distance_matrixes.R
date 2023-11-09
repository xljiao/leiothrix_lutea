library(geosphere)
library(Cairo)
library(sp)
library(raster)

# calculate the distance of two sites using geographical coordinates

location<-read.csv("location.csv",header = TRUE)
colnames(location)[1]<-"ID"
muer.lonlat = cbind(location$log, location$lat)
muer.dists = distm(muer.lonlat, fun=distVincentyEllipsoid);
rownames(muer.dists) = location$ID;
colnames(muer.dists) = location$ID;
write.table(muer.dists, 'matrix_distance.csv', quote=FALSE, sep=",");

# extract bioclimatic variables of each sites
LonLatData <- location[,c(2,3)]
#input bioclimatic variables
files <- list.files("current_bio_30s",pattern='asc',full.names=TRUE)
Grids <- raster::stack(files)
ID<-location[,c(1)]
VariablesAtlocations <- raster::extract(Grids,LonLatData)
finalfile<-as.data.frame(cbind(ID,LonLatData,VariablesAtlocations))
colnames(finalfile) <- c("species","longitude","latitude",
                         colnames(VariablesAtlocations))
write.table(finalfile, 'finalfile.csv', quote=FALSE, sep=",")

#PCA analysis of sampling sites
location_bio<-read.csv("finalfile.csv",header = TRUE)
location_bio <- location_bio[,-c(1,2,3)]
pca1<-princomp(location_bio,cor = TRUE)
plot(pca1)
biplot(pca1)
pca1$loadings
summary(pca1, loadings = TRUE)
pca_loadings <- pca1$loadings
write.table(pca_loadings, 'pca_loadings.csv', quote=FALSE, sep=",")
pca1$scores[,1:2]
pca1_scores <- pca1$scores[,1]
write.table(pca1_scores, 'pca1_scores.csv', quote=FALSE, sep=",")
pca2_scores <- pca1$scores[,2]
write.table(pca2_scores, 'pca2_scores.csv', quote=FALSE, sep=",")

### calculate the euclidean distance of PC1 and PC2
pc1 <- read.csv(file="pca1_scores.csv",header=T,row.names=1)
edist_pc1 <- as.matrix(dist(pc1))
write.table(edist_pc1, 'matrix_bio_pc1.csv', quote=FALSE, sep=",")

pc2 <- read.csv(file="pca2_scores.csv",header=T,row.names=1)
edist_pc2 <- as.matrix(dist(pc2))
write.table(edist_pc2, 'matrix_bio_pc2.csv', quote=FALSE, sep=",")


