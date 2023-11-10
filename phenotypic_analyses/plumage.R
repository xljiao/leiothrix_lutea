library(pavo)

pluamge_specs <- getspec("pluamge",ext = "txt",lim = c(300,700),decimal = ".",
                       subdir = TRUE, subdir.names = FALSE)

# Calculate the average of 5 measurements for each part
pluamge_mspecs <- aggspec(pluamge_specs,by = 5)

# Use 0.2 as a threshold to eliminate noise
plotsmooth(crown_specs,minsmooth = 0.05,maxsmooth = 0.5,curves = 4,ask = FALSE) 
pluamge.sm <- procspec(crown_mspecs,opt="smooth",span = 0.2) 

# We used the "vismodel" function to convert spectral data into quantum catches of the bird retina under the tetrachromatic visual system, using the blue tit as a visual model 
all_specs <- vismodel(pluamge.sm,visual = c("bluetit"),relative = FALSE)
write.table(all_specs,file = "all_specs",sep=",") 
