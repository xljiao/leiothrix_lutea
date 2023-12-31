### MLPE
library(ResistanceGA)
library(raster)
library(sp)
library(lme4)
fst <- as.matrix(read.csv("matrix_fst2.csv"))
res <- as.matrix(read.csv("matrix_resistance.csv"))
geo <- as.matrix(read.csv("matrix_distance.csv"))
pc1 <- as.matrix(read.csv("matrix_bio_pc1.csv"))
pc2 <- as.matrix(read.csv("matrix_bio_pc2.csv"))
MLPE.lmm(res, fst,REML = FALSE, ID = NULL, ZZ = NULL, scale = TRUE)
MLPE.lmm(geo, fst,REML = FALSE, ID = NULL, ZZ = NULL, scale = TRUE)
MLPE.lmm(pc1, fst,REML = FALSE, ID = NULL, ZZ = NULL, scale = TRUE)
MLPE.lmm(pc2, fst,REML = FALSE, ID = NULL, ZZ = NULL, scale = TRUE)

## MMRR
library(sp)
library(raster)
library(ecodist)

alldata <- read.csv("all.csv")
MRM(fst2 ~ geo+pc1+pc2+res, data=alldata, nperm=10000,method = "linear",mrank = TRUE)
MRM(fst2 ~ geo, data=alldata, nperm=10000,method = "linear",mrank = TRUE)
MRM(fst2 ~ pc1, data=alldata, nperm=10000,method = "linear",mrank = TRUE)
MRM(fst2 ~ pc2, data=alldata, nperm=10000,method = "linear",mrank = TRUE)
MRM(fst2 ~ res, data=alldata, nperm=10000,method = "linear",mrank = TRUE)
