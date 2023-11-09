## create sampling sites raster, the grid size needs to be the same as the resistance surface grid size
testa <- raster("leiothrixlutea_avg.asc")  # resistance surface raster layer
k <- read.csv("location.csv") #location
cor <- k[,c(2,3)] #extract the geographical coordinates of sampling sites
num <- cellFromXY(testa,cor)  #Find the grid number of the coordinate point
a <- as.data.frame(testa)  # Assign the drag surface raster file to a
a$leiothrixlutea_avg <- c(rep(-9999,nrow(a))) # test is the column name. The name is based on the file name. The test column (only this column) is changed to -9999
a[num,1] <- k$value  #Assigns the grid number of the coordinate point to the number set for the site-point in the coordinate point file
values(testa) <- a$leiothrixlutea_avg #Assigns the modified raster file value to the original raster file
plot(testa) # can look directly at the current value distribution of the grid
writeRaster(testa,"location_newn.asc") #Export to raster file
