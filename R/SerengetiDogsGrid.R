#
## Get information about human & dog populations in cells making up Mara grid
#____________________________

rm(list=ls())

library(rgdal)
library(rgeos)
library(RColorBrewer)
library(raster)

options(stringsAsFactors=F) 



## Load data files
#--------------

startYear <- 2000
startYear_analysis <- 2002
endYear <- 2040
res<-1

##Load shapefile
SD_vill <- readOGR("output/SD_vill","SD_vill") #(version created in SerengetDogsVillage.R)

## Read in village populations
humanMatVill <- as.matrix(read.csv(paste("output/humanPopulationByVillageMonth_Jan",startYear,"_Dec",endYear,".csv",sep=""),row.names = 1,header=F))
dogMatVill <- as.matrix(read.csv(paste("output/dogPopulationByVillageMonth_Jan",startYear,"_Dec",endYear,".csv",sep=""),row.names = 1,header=F))

# Rasters
villGrid <- raster(paste("output/villGrid_",res,"km.grd",sep=""))
distGrid <- raster(paste("output/cellGrid_",res,"km.grd",sep=""))





## Create dog density grid from census data
#--------------

# Code for this step is commented out here for reference, but original census
# data file with coordinates of households is not included in this public repo.
# The grids resulting from this step are included and read in as a raster at the
# end of this code section.

# ## Load Serengeti census data and convert to spatial points dataframe
# Census <- read.csv("data/SDcompiled.csv")
# Census <- SpatialPointsDataFrame(Census[,c(10,9)],Census, proj4string = CRS("+proj=longlat"))
# Census <- spTransform(Census, SD_vill@proj4string)
# 
# ## Find grid cell associated with each household
# Census$cellID <- raster::extract(distGrid,Census)
# 
# ##Assign any census points that landed in NA cells to the closest non-NA cell
# hhToAssign <- which(is.na(Census$cellID))
# Census$cellID[hhToAssign] <- apply(X = matrix(Census@coords[hhToAssign,],ncol=2), MARGIN = 1, 
#                        FUN = function(xy) distGrid[which.min(replace(distanceFromPoints(distGrid, xy), is.na(distGrid), NA))])
# 
# ## Human population density grid 
# humans<-rowsum(Census$humansTotal,Census$cellID)
# humanGrid <- distGrid
# humanGrid[which(!is.na(humanGrid[]))]<-0
# humanGrid[which(!is.na(humanGrid[]))[as.numeric(rownames(humans))]] <- humans
# 
# ## Dog population density grid 
# dogs<-rowsum(Census$dogsTotal,Census$cellID)
# dogGrid <- distGrid
# dogGrid[which(!is.na(dogGrid[]))]<-0
# dogGrid[which(!is.na(dogGrid[]))[as.numeric(rownames(dogs))]] <- dogs
# writeRaster(dogGrid,file=paste("data/GIS/dogGrid_1km.grd",sep=""),overwrite=T)

# Read in final dog and human grids
dogGrid <- raster(paste("data/GIS/dogGrid_1km.grd",sep=""))
humanGrid <- raster(paste("data/GIS/humanGrid_1km.grd",sep=""))

## Number of cells in each village that are occupied
for(i in 1:nrow(SD_vill)){
  SD_vill$occ_Dog[i] <- length(which(dogGrid[which(villGrid[]==i)]>0))
  SD_vill$cells_occ[i] <- length(which(humanGrid[which(villGrid[]==i)]>0))
}
c(rowsum(dogGrid[which(!is.na(villGrid[]))],villGrid[which(!is.na(villGrid[]))])/SD_vill$cells_occupied)
writeOGR(SD_vill, dsn="output/SD_vill", "SD_vill", driver="ESRI Shapefile",overwrite_layer=T)



## Get dog and human population by cell and month
#--------------

dogMat <- matrix(NA,ncol=ncol(dogMatVill),nrow=length(dogGrid[which(!is.na(dogGrid[]))]))
humanMat <- matrix(NA,ncol=ncol(dogMatVill),nrow=length(dogGrid[which(!is.na(dogGrid[]))]))

for(i in 1:nrow(dogMatVill)){
  
  cells <- which(villGrid[]==i)
  probs_dog <- dogGrid[cells]/sum(dogGrid[cells])
  probs_human <- humanGrid[cells]/sum(humanGrid[cells])
  
  cells <- which(villGrid[which(!is.na(villGrid[]))]==i)
  dogMat[cells,] <- rep(dogMatVill[i,],each=length(cells)) * probs_dog
  humanMat[cells,] <- rep(humanMatVill[i,],each=length(cells)) * probs_human

}


## Save
write.table(humanMat,file=paste("output/humanPopulationByCellMonth_Jan",startYear,"_Dec",endYear,"_",res,"kmGrid.csv",sep=""), sep=",", row.names = F, col.names = F)
write.table(dogMat,file=paste("output/dogPopulationByCellMonth_Jan",startYear,"_Dec",endYear,"_",res,"kmGrid.csv",sep=""), sep=",", row.names = F, col.names = F)
write.table(humanMat[,((startYear_analysis-startYear)*12+1):ncol(humanMat)],file=paste("output/humanPopulationByCellMonth_Jan",startYear_analysis,"_Dec",endYear,"_",res,"kmGrid.csv",sep=""), sep=",", row.names = F, col.names = F)
write.table(dogMat[,((startYear_analysis-startYear)*12+1):ncol(dogMat)],file=paste("output/dogPopulationByCellMonth_Jan",startYear_analysis,"_Dec",endYear,"_",res,"kmGrid.csv",sep=""), sep=",", row.names = F, col.names = F)




## Plot
#--------------

colours <- colorRampPalette(c("white",brewer.pal(8,"YlOrRd")))(100)

par(mfrow=c(1,2),mar=c(0,0,0,6))
popGrid<-distGrid
popGrid[which(!is.na(villGrid[]))]<-humanMat[,1]/(res^2)
popGrid <- log10(popGrid)
popGrid[which(popGrid[]<0)]<-0
plot(SD_vill)
plot(popGrid,add=T,col=colours,breaks=seq(0,max((popGrid[]),na.rm=T),length.out=100),legend=F)
plot(SD_vill,add=T)
plot(popGrid, breaks=seq(0,max((popGrid[]),na.rm=T),length.out=100),
     legend.only=T, add=T,col=colours,
     legend.args=list(text=bquote("Humans/"~.(res^2) * "km"^2), side=4, font=2, line=3, cex=1.2),
     axis.args=list(at=c(log10(c(1,10,100,1000))),labels=c(expression(""<=1),"10","100","1000")),cex.axis=0.8,
     smallplot=c(0.75,0.76, .25,.75))

popGrid<-distGrid
popGrid[which(!is.na(villGrid[]))]<-dogMat[,1]/(res^2)
popGrid <- log10(popGrid)
popGrid[which(popGrid[]<0)]<-0
plot(SD_vill)
plot(popGrid,add=T,col=colours,breaks=seq(0,max((popGrid[]),na.rm=T),length.out=100),legend=F)
plot(SD_vill,add=T)
plot(popGrid, breaks=seq(0,max((popGrid[]),na.rm=T),length.out=100),
     legend.only=T, add=T,col=colours,
     legend.args=list(text=bquote("Dogs/"~.(res^2) * "km"^2), side=4, font=2, line=3, cex=1.2),
     axis.args=list(at=c(log10(c(1,10,100,1000))),labels=c(expression(""<=1),"10","100","1000")),cex.axis=0.8,
     smallplot=c(0.75,0.76, .25,.75))





