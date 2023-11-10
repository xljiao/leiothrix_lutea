library(rEEMSplots)
library(rgdal)
library(rworldmap)
library(rworldxtra)

PROJ_LIB=~/EEMS/share/proj

projection_none <- "+proj=longlat +datum=WGS84"

## 所有10次合并

eems.plots(mcmcpath = c("eems_chain1","eems_chain2","eems_chain3","eems_chain4","eems_chain5","eems_chain6","eems_chain7","eems_chain8","eems_chain9","eems_chain10"),
	plotpath = paste0("~/plot/", "-chain1-10"),
	longlat = TRUE,
	projection.in = projection_none,
	projection.out = projection_none,
	add.demes =TRUE,
	add.map = TRUE,
	out.png = FALSE,
	lwd.map = 1)
