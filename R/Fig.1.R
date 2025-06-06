rm(list=ls())

library(raster)
library(viridis)
library(RColorBrewer)
library(rgeos)
library(prettymapr)
library(rgdal)


## Load data
#______________________

# Serengeti shapefile
SD_vill <- readOGR("Output/SD_vill","SD_vill") 
SD_vill<-SD_vill[order(SD_vill$Vill_2012),]
SD_outline <- gUnaryUnion(SD_vill)# get district outline
SD_outline<-gBuffer(SD_outline,width=1) # get rid of a few tiny holes

# Protected areas
PAs <- readOGR("data/GIS/ProtectedAreas","TZprotected_areas", p4s = ("+proj=utm +zone=37 +south +ellps=clrk80 +towgs84=-160,-6,-302,0,0,0,0 +units=m +no_defs"))
PAs <- spTransform(PAs,SD_outline@proj4string)

# Tanzania outline
TzOutline <- readOGR("data/GIS/TZ_Outline_2012","TZ_Outline_2012")
TzOutline <- spTransform(TzOutline,SD_vill@proj4string)

# Serengeti grid
res<-1
villGrid <- raster(paste("Output/villGrid_",res,"km.grd",sep=""))

# Dog population
years <- 2002:2022
dogs <- as.matrix(read.csv("Output/dogPopulationByCellMonth_Jan2002_Dec2040_1kmGrid.csv",header=F))
dogs_vill <- as.matrix(read.csv("Output/dogPopulationByVillageMonth_Jan2002_Dec2022.csv",header=F,row.names = 1))
dog_dens_vill <- dogs_vill/SD_vill$cells_occ
dogs <- dogs[,1:(12*length(years))]

# Human population
humans <- as.matrix(read.csv("Output/humanPopulationByCellMonth_Jan2002_Dec2040_1kmGrid.csv",header=F))
humans <- humans[,1:(12*length(years))]



## Plot
#______________________

pdf("Figs/Fig1.pdf",width=7,height=4.5)

colours <- colorRampPalette(c("white",brewer.pal(8,"YlOrRd")))(100)
cex.axis <- 0.7
cex.lab <- 0.8

# Tanzania inset
par(fig=c(0,0.20,0.8,1),mar=c(0,0,0,0),bg="transparent")
plot(TzOutline)
rect(121654.2-30000,8697442.6-30000,1317214+30000,9891260+30000)
text(TzOutline,"Tanzania",cex=cex.lab)
plot(SD_outline,add=T,col="blue3",border=NA)
plot(PAs[which(PAs$SP_ID=="2"),],add=T,col="grey",border=NA)
legend(-650000,9891260,legend="A",text.font = 2,bty="n",xpd=F)

# Study area map with dog population
par(fig=c(0,0.5,0.45,1),mar=c(0,0,0,0),new=T)
popGrid<-villGrid
popGrid[which(!is.na(villGrid[]))]<-dogs[,ncol(dogs)]
popGrid <- log10(popGrid)
popGrid[which(popGrid[]<0)]<-0
plot(SD_outline,xlim=c(637188-30000,707444.2),ylim=c(9754421-8000,9837984.4))
plot(popGrid,add=T,col=colours,breaks=seq(0,max((popGrid[]),na.rm=T),length.out=100),legend=F)
par(fig=c(0,0.5,0.45,1),mar=c(0,0,0,0),new=T)
plot(SD_vill,add=T,border="grey60",lwd=0.5)
plot(PAs[which(PAs$SP_ID=="2"),],add=T,col="grey",border=NA)
plot(SD_outline,add=T)
addscalebar(plotunit="m",plotepsg = 32736,htin = 0.05)
legend(640000,9754421,"Serengeti National Park",bty="n",cex=cex.lab)
plot(popGrid, breaks=seq(0,max((popGrid[]),na.rm=T),length.out=100),
     legend.only=T, add=T,col=colours,
     legend.args=list(text=bquote("Dogs/"~.(res^2) * "km"^2), side=4, line=1.7, cex=cex.lab),
     axis.args=list(at=c(log10(c(1,10,100,1000))),labels=c(expression(""<=1),"10","100","1000"),hadj=0.5,cex.axis=cex.axis),
     smallplot=c(0.9,0.92, .25,.75))

# addnortharrow("topright",scale = 0.8)
par(fig=c(0.6,1,0.45,1),mar=c(3,3,1.5,1),new=T)
plot(colSums(dogs)~c(1:ncol(dogs)),type="l",col="navy",axes=F,lwd=2,ylab="",xlab="",ylim=c(0,90000))
box(bty="l")
axis(1,at=seq(1,length(2002:2022)*12,48),labels=paste(seq(2002,2022,4)),cex.axis=cex.axis,padj=-1)
axis(2,cex.axis=cex.axis,padj=1.2,at=c(0,40000,80000),labels=format(c(0,40000,80000),scientific = F,big.mark=",",trim=T))
mtext("Estimated dog population",side=2,line=1.6,cex=cex.lab)
mtext("Date",side=1,line=1.5,cex=cex.lab)
legend("topleft",legend="B",text.font = 2,bty="n")

# Blank out bottom portion
par(fig=c(0,1,0,0.45),mar=c(0,0,0,0),new=T)
plot(SD_outline,bg="white",border="white")

# Human dog ratio
par(fig=c(0,0.5,0,0.42),mar=c(0,0,0,4),new=T)
breaks=seq(min(SD_vill$HDR),max(SD_vill$HDR),length.out=100)
plot(SD_vill,border="grey60",lwd=0.5, col=colours[findInterval((SD_vill$HDR),breaks,all.inside=T)])
plot(SD_outline,add=T)
legend("topleft",legend="C",text.font = 2,bty="n",xpd=F)
grid <- raster(extent(SD_vill),crs=SD_vill@proj4string);res(grid) <- 1000;grid[]<-5
plot(grid, 
     breaks=breaks,legend.only=T,col=colours,
     legend.args=list(text="Human:dog ratio", side=4, line=1.5, cex=cex.lab),
     axis.args=list(at=pretty(c(min(SD_vill$HDR),max(SD_vill$HDR)),n = 4),cex.axis=cex.axis),
     smallplot=c(0.70,0.72, .2,.8))

# Village density
par(fig=c(0.5,1,0,0.42),mar=c(0,0,0,4),new=T)
breaks=seq(0,max((dog_dens_vill[]),na.rm=T),length.out=100)
plot(SD_vill,border="grey60",lwd=0.5, col=colours[findInterval((dog_dens_vill[,ncol(dogs)]),breaks,all.inside=T)])
plot(SD_outline,add=T)
legend("topleft",legend="D",text.font = 2,bty="n",xpd=F)
grid <- raster(extent(SD_vill),crs=SD_vill@proj4string);res(grid) <- 1000;grid[]<-1
plot(grid, 
     breaks=breaks,legend.only=T,col=colours,
     legend.args=list(text=bquote("Dogs/"~.(res^2) * "km"^2), side=4, line=2.2, cex=cex.lab),
     axis.args=list(at=seq(0,250,50),cex.axis=cex.axis),
     smallplot=c(0.70,0.72, .2,.8))

dev.off()


pdf("Figs/Fig1_unlogged.pdf",width=7,height=2.5)

colours <- colorRampPalette(c("white",brewer.pal(8,"YlOrRd")))(100)
cex.axis <- 0.7
cex.lab <- 0.8

# Tanzania inset
par(fig=c(0,0.20,0.6,1),mar=c(0,0,0,0),bg="transparent")
plot(TzOutline)
text(TzOutline,"Tanzania",cex=cex.lab)
plot(SD_outline,add=T,col="blue3",border=NA)
plot(PAs[which(PAs$SP_ID=="2"),],add=T,col="grey",border=NA)
legend(-450000,9891260,legend="A",text.font = 2,bty="n",xpd=F)

# Study area map with dog population
par(fig=c(0,0.5,0,1),mar=c(0,0,0,0),new=T)
popGrid<-villGrid
popGrid[which(!is.na(villGrid[]))]<-dogs[,ncol(dogs)]
plot(SD_outline,xlim=c(637188-30000,707444.2),ylim=c(9754421-8000,9837984.4))
plot(popGrid,add=T,col=colours,breaks=seq(0,max((popGrid[]),na.rm=T),length.out=100),legend=F)
par(fig=c(0,0.5,0,1),mar=c(0,0,0,0),new=T)
plot(SD_vill,add=T,border="grey60",lwd=0.5)
plot(PAs[which(PAs$SP_ID=="2"),],add=T,col="grey",border=NA)
plot(SD_outline,add=T)
addscalebar(plotunit="m",plotepsg = 32736,htin = 0.05)
legend(640000,9754421,"Serengeti National Park",bty="n",cex=cex.lab)
plot(popGrid, breaks=seq(0,max((popGrid[]),na.rm=T),length.out=100),
     legend.only=T, add=T,col=colours,
     legend.args=list(text=bquote("Dogs/"~.(res^2) * "km"^2), side=4, line=1.7, cex=cex.lab),
     axis.args=list(at=pretty(c(0,max(popGrid[],na.rm=T))),hadj=0.5,cex.axis=cex.axis),
     smallplot=c(0.9,0.92, .25,.75))

# addnortharrow("topright",scale = 0.8)
par(fig=c(0.6,1,0,1),mar=c(3,3,1.5,1),new=T)
plot(colSums(dogs)~c(1:ncol(dogs)),type="l",col="navy",axes=F,lwd=2,ylab="",xlab="",ylim=c(0,90000))
box(bty="l")
axis(1,at=seq(1,length(2002:2022)*12,48),labels=paste(seq(2002,2022,4)),cex.axis=cex.axis,padj=-1)
axis(2,cex.axis=cex.axis,padj=1.2,at=c(0,40000,80000),labels=format(c(0,40000,80000),scientific = F,big.mark=",",trim=T))
mtext("Estimated dog population",side=2,line=1.6,cex=cex.lab)
mtext("Date",side=1,line=1.5,cex=cex.lab)
legend("topleft",legend="B",text.font = 2,bty="n")


dev.off()



# Some summaries for MS
#______________________

range(popGrid[],na.rm=T)

n <- length(popGrid[which(!is.na(popGrid[]))])
n2 <- length(popGrid[which(popGrid[]>0)])
round((length(popGrid[which(popGrid[]!=0)])/length(which(!is.na(popGrid[]))))*100)

mean(dogs[,ncol(dogs)])
sum(dogs[,ncol(dogs)])/n2

100*length(which(popGrid[]>100))/n
100*length(which(popGrid[]>100))/n2
100*sum(popGrid[which(popGrid[]>100)],na.rm=T)/sum(popGrid[],na.rm=T)

100*length(which(popGrid[]>200))/n
100*length(which(popGrid[]>200))/n2

gArea(SD_vill)/1e6

range(colSums(dogs))

round(range(dog_dens_vill[,ncol(dogs)]))
range(SD_vill$HDR)

as.numeric(colSums(humans)/colSums(dogs))
