rm(list=ls())

library(raster)
library(rgdal)

options(stringsAsFactors=F) 
res <- 1 # cells are resxres km in size


## Load in shapefiles
SD_vill <- readOGR("data/GIS","SD_Villages_2012_From_HHS_250m_Lattice_UTM") 
SD_vill<-SD_vill[order(SD_vill$Vill_2012),]


## Create grid
grid <- raster(extent(SD_vill),crs=SD_vill@proj4string)
res(grid) <- res*1000
gridPoints <- SpatialPoints(rasterToPoints(grid), proj4string = SD_vill@proj4string)


## Create 2012 village grid
SD_vill$VillID <- 1:nrow(SD_vill)
values <- over(gridPoints,SD_vill)$VillID
villGrid <- grid
villGrid[] <- values
plot(villGrid)

## Create 2012 ward grid
SD_vill$wardID <- match(SD_vill$Ward_2012,sort(unique(SD_vill$Ward_2012)))
values <- over(gridPoints,SD_vill)$wardID
wardGrid <- grid
wardGrid[] <- values
plot(wardGrid)

## District grid
cellGrid <- villGrid
cellGrid[which(!is.na(villGrid[]))] <- 1:length(which(!is.na(villGrid[])))
plot(cellGrid)

##Save grids
writeRaster(villGrid,file=paste("output/villGrid_",res,"km.grd",sep=""),overwrite=T)
writeRaster(cellGrid,file=paste("output/cellGrid_",res,"km.grd",sep=""),overwrite=T)
write.table(as.matrix(cellGrid),paste("output/cell_matrix_",res^2,"kmsq.csv",sep=""),row.names=F,col.names=F,sep=",")

## Data frame with info on each cell
cellData<-data.frame(cellID=1:length(which(!is.na(cellGrid[]))))
cellData$Region <- "Mara"
cellData$District <- "Serengeti"
cellData$Ward <- sort(unique(SD_vill$Ward_2012))[(wardGrid[which(!is.na(wardGrid@data@values))])]
cellData$WardID <- wardGrid[which(!is.na(wardGrid@data@values))]
cellData$Village <- SD_vill$Vill_2012[(villGrid[which(!is.na(villGrid@data@values))])]
cellData$VillageID <- villGrid[which(!is.na(villGrid@data@values))]
write.table(cellData,paste("output/SerengetiCellData_",res^2,"kmsq.csv",sep=""),row.names=F,sep=",")



