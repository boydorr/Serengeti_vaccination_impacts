
rm(list=ls())

library(raster)
library(viridis)
library(RColorBrewer)
library(rgeos)
library(animation)
library(rgeos)
library(rgdal)


## Load data
#______________________

years <- 2002:2022

# Serengeti village shapefile
SD_vill <- readOGR("output/SD_vill","SD_vill") 
SD_vill<-SD_vill[order(SD_vill$Vill_2012),]
SD_outline <- gUnaryUnion(SD_vill)# get district outline
SD_outline<-gBuffer(SD_outline,width=1) # get rid of a few tiny holes
mean(sqrt(gArea(SD_vill,byid = T)))/1000
mean(sqrt(SD_vill$cells_occ))

# Coverage estimates
vax_vill <- as.matrix(read.csv("output/vaccinationCoverageByVillageMonth_Jan2002_Dec2022.csv",header = F,row.names = 1)) # village level
vax_dist <- as.matrix(read.csv("output/districtVaccinationCoverage_Jan2002_Dec2022.csv",header = T))[,1] # district level
vax_vill_annual <- matrix(NA,nrow=nrow(vax_vill),ncol=length(years))
vax_dist_annual <- rep(NA,length(years))
for(i in 1:length(years)){
  vax_vill_annual[,i] <- rowMeans(vax_vill[,(1:12) + (i-1)*12])
  vax_dist_annual[i] <- mean(vax_dist[(1:12) + (i-1)*12])
}
names(vax_dist_annual) <- years

# Dog numbers 
dogs <- as.matrix(read.csv("output/dogPopulationByVillageMonth_Jan2002_Dec2040.csv",row.names = 1,header=F))
dogs <- dogs[,1:(length(years)*12)]
dogs_dist <- colSums(dogs)
dogs_dist_annual <- rep(NA,length(years))
for(i in 1:length(years)){
  dogs_dist_annual[i] <- round(mean(dogs_dist[(1:12) + (i-1)*12]))
}

# Human numbers 
humans <- as.matrix(read.csv("output/humanPopulationByVillageMonth_Jan2002_Dec2040.csv",row.names = 1,header=F))
humans <- humans[,1:(length(years)*12)]
humans_dist <- colSums(humans)
humans_dist_annual <- rep(NA,length(years))
for(i in 1:length(years)){
  humans_dist_annual[i] <- round(mean(humans_dist[(1:12) + (i-1)*12]))
}

# campaign coverage and numbers vaccinated
n_dogs_vax <- read.csv("output/dogsVaccinatedByVillageYear.csv")
n_dogs_vax_dist <- colSums(n_dogs_vax)
percent_dogs_vax <- c(round(as.matrix(read.csv("output/vcDistByYear.csv",header = F))*100)) # annual district campaign coverages
villages_vax <- read.csv("output/VillagesVaccinatedByYear.csv")
n_villages_vax <- colSums(villages_vax)
percent_villages_vax <- round((n_villages_vax/nrow(SD_vill))*100)
campaign_cov_vill <- as.matrix(read.csv("output/vcVillByYear.csv",header=F))

# Dog case numbers
dog_cases_dist <- as.matrix(read.csv("output/Serengeti_monthly_rabid_dogs_2002-01-01_2022-12-31.csv",header=F))
dog_cases_vill <- as.matrix(read.csv("output/Serengeti_monthly_rabid_dogs_village_2002-01-01_2022-12-31.csv",header=F,row.names = 1))
dog_cases_CT <- read.csv("output/Serengeti_rabid_dogs.csv")
dog_cases_CT$year <- (dog_cases_CT$month-1)%/%12+2002

# Human bites by month
monthlyBites <- read.csv("output/Serengeti_monthly_human_exposures_2002-01-01_2022-12-31.csv",header=F)$V1 # exposures by all species
exposures_CT <- read.csv("output/Serengeti_exposed_humans.csv")
exposures_CT$year <- (exposures_CT$month-1)%/%12+2002

# Human bites by dogs by month
monthlyBitesDogs <- read.csv("output/Serengeti_monthly_human_exposures_dogs_2002-01-01_2022-12-31.csv",header=F)$V1

# Human deaths by month
monthlyDeaths <- read.csv("output/Serengeti_monthly_human_deaths_2002-01-01_2022-12-31.csv",header=F)$V1 
deaths_CT <- read.csv("output/Serengeti_human_deaths.csv")
deaths_CT$year <- (deaths_CT$month-1)%/%12+2002



## Manuscript summaries
#______________________

# Heterogeneity in coverage
var_cov <- rep(NA,length(vax_dist))
for(i in 1:length(vax_dist)){var_cov[i] <- sum((dogs[,i]/sum(dogs[,i]))*(vax_vill[,i]-vax_dist[i])^2)}
sd_cov <- sqrt(var_cov)
sd_cov_annual <- rep(NA,length(vax_dist_annual))
for(i in 1:length(years)){
  sd_cov_annual[i] <- mean(sd_cov[(1:12) + (i-1)*12])
}
range(sd_cov)
names(sd_cov_annual) <- 2002:2022
sort(round(sd_cov_annual,2))
which.max(sd_cov_annual)

# Some campaign coverage summaries for MS
n_dogs_vax_dist
range(n_dogs_vax_dist[-1])
range(percent_dogs_vax[-1])
mean(percent_dogs_vax[-1])
range(campaign_cov_vill[,-1])*100
range(campaign_cov_vill[which(campaign_cov_vill!=0)])*100
mean(campaign_cov_vill[,-1])*100
mean(campaign_cov_vill[which(campaign_cov_vill!=0)])*100
100*length(which(campaign_cov_vill[,-1]>=0.7))/length(c(campaign_cov_vill[,-1]))
100*length(which(campaign_cov_vill[,-1]>=0.7))/length(c(campaign_cov_vill[which(campaign_cov_vill!=0)]))

# Some coverage summaries for MS
range(vax_dist)*100
mean(vax_dist)*100
range(vax_dist[(12*3+1):length(vax_dist)])*100
diff(range(vax_dist[(12*3+1):length(vax_dist)])*100)
diff(range(vax_dist_annual[3:length(vax_dist_annual)]))
quantile(range(vax_dist_annual[3:length(vax_dist_annual)]))
vax_dist[1]*100
range(vax_vill)*100
range(vax_vill[which(vax_vill!=0)])*100
length(which(vax_dist<0.2))/length(vax_dist) # how often does coverage fall below critical point
for(y in 2002:2022){
  print(c(y,100*length(which(vax_dist[((y-2002)*12+1):length(vax_dist)]<0.2))/length(vax_dist[((y-2002)*12+1):length(vax_dist)])))
}

# dog case summaries for manuscript
nrow(dog_cases_CT)
table(dog_cases_CT$year)
range(table(dog_cases_CT$year))
(table(dog_cases_CT$year)/dogs_dist_annual)*100
range((table(dog_cases_CT$year)/dogs_dist_annual)*100)
annual_incidence <- (table(dog_cases_CT$year)/(dogs_dist_annual/1000))
annual_incidence
annual_incidence[which(annual_incidence>4)]
range((table(dog_cases_CT$year)/(dogs_dist_annual/1000)))
range(dog_cases_dist)
range(dog_cases_vill)
range(rowSums(dog_cases_vill))
length(which(rowSums(dog_cases_vill)==0))
SD_vill$Vill_2012[which(rowSums(dog_cases_vill)==0)]
plot(SD_vill)
plot(SD_vill[which(rowSums(dog_cases_vill)==0),],add=T,col=2)
plot(SD_vill[which(rowSums(dog_cases_vill)==1),],add=T,col=3)
plot(SD_vill[which(rowSums(dog_cases_vill)==2),],add=T,col=4)
hist(rowSums(dog_cases_vill),breaks=seq(-0.5,186+0.5,1))
head(sort(rowSums(dog_cases_vill)),10)
tail(sort(rowSums(dog_cases_vill)),10)

# exposure summaries for manuscript
nrow(exposures_CT)
table(exposures_CT$year)
range(table(exposures_CT$year))
(table(exposures_CT$year)/humans_dist_annual)*1000
range(monthlyBites)
which.max(monthlyBites)

# death  summaries for MS
nrow(deaths_CT)
table(deaths_CT$year)
deaths_annual <- rep(0,length(years))
names(deaths_annual) <- years
deaths_annual[names(table(deaths_CT$year))] <- table(deaths_CT$year)
(deaths_annual/humans_dist_annual)*100000
max(table(deaths_CT$year))
c(years)[which(!c(years)%in%names(table(deaths_CT$year)))]
table(deaths_CT$PEP.1)
table(deaths_CT$PEP.2)
table(deaths_CT$PEP.1,deaths_CT$PEP.2)



# Comparison of bites by dogs and all bites
#______________________

sum(monthlyBitesDogs)/sum(monthlyBites)
hist(monthlyBites-monthlyBitesDogs,breaks=seq(-0.5,10.5,1))
length(which(exposures_CT$Attacking.species=="Domestic dog"))/nrow(exposures_CT)
exposures_CT$species <- exposures_CT$Attacking.species
exposures_CT$species[which(exposures_CT$species=="Other")] <- exposures_CT$Other.attacking.species[which(exposures_CT$species=="Other")]
# round((table(exposures_CT$species)/nrow(exposures_CT))*100,1)
# round((table(exposures_CT$species[which(exposures_CT$species!="Domestic dog")])/nrow(exposures_CT[which(exposures_CT$species!="Domestic dog"),]))*100,1)
not_dog <- which(!exposures_CT$species%in%c("Domestic dog","","dog/wildlife?"))
n_not_dog <- length(not_dog)
sort(round((table(exposures_CT$species[not_dog])/n_not_dog)*100,1))
length(which(exposures_CT$Attacking.species=="Wildlife: Jackal"))/n_not_dog
round((sum(grepl("Livestock",exposures_CT$species))/n_not_dog)*100,1)
round((sum((!(grepl("Livestock",exposures_CT$species)|exposures_CT$species=="Human"|exposures_CT$species%in%c("Domestic dog","","dog/wildlife?","Cat","Wildlife: Jackal"))))/n_not_dog)*100,1)

pdf("Figs/ExposuresByDogsVsOtherSpp.pdf",width=7, height=4)
cex.axis <- 0.7
cex.lab <- 0.8
par(mar=c(3,2.5,1,0.5),mfrow=c(1,1))
plot(monthlyBitesDogs,axes=F,col="orange",ylim=c(0,max(monthlyBitesDogs)),type="l",lwd=2,ylab="",xlab="")
lines(monthlyBites-monthlyBitesDogs,axes=F,col="navy",type="l",lwd=2)
axis(2,at=seq(0,30,10),cex.axis=cex.axis,padj=1)
axis(1,at=seq(1,length(years)*12,24),labels=paste(seq(min(years),max(years),2)),cex.axis=cex.axis,padj=-1)
mtext("Human exposures",line=1.5,side=2,cex=cex.lab)
mtext("Date",side=1,cex=cex.lab,line=1.5)
legend("topright",c("Exposures by dogs","Exposures by all other species"),col=c("orange","navy"),lwd=2,lty=1,bty="n",cex=cex.lab)
box(bty="l")
dev.off()



# Supplementary table 1
#______________________

ST1 <- rbind(n_dogs_vax_dist, percent_dogs_vax, percent_villages_vax)
colnames(ST1) <- years
write.csv(ST1, file="output/vaccination_table.csv")



## Fig. 2
#______________________

pdf("Figs/Fig.2.pdf",width=7, height=8.75)
cex.axis <- 0.7
cex.lab <- 0.8
cex.pt <- 0.2



## Campaign completeness
#------------

par(fig=c(0,0.5,0.79,1))
par(mar=c(2,3.5,0.9,0))

percent_villages_vax
b<-barplot(percent_villages_vax,names=NA,axes=F,space=0.2,col="goldenrod1",border=NA,cex.names=cex.axis,ylim=c(0,100))
axis(2,cex.axis=cex.axis,padj=1.2)#,at=seq(1,length(years)*12,24),labels=paste(seq(min(years),max(years),2)))
axis(2,cex.axis=cex.axis,at=0,padj=1.2)#,at=seq(1,length(years)*12,24),labels=paste(seq(min(years),max(years),2)))
axis(1,at=b[c(T,F,F,F)],labels=paste(seq(min(years),max(years),4)),cex.axis=cex.axis,padj=-2)
mtext("Campaign completeness\n(% villages vaccinated)",side=2,line=1.6,cex=cex.lab)
mtext("Date",side=1,line=1,cex=cex.lab)
box(bty="l")
legend(-3,100,legend="A",text.font = 2,bty="n")



## Vaccination numbers
#------------

par(fig=c(0.5,1,0.79,1),new=T)
par(mar=c(2,3,0.9,0.5))

b<-barplot(dogs_dist_annual,names=NA,axes=F,space=0.2,col="navy",border=NA,cex.names=cex.axis,ylim=c(0,max(dogs_dist)))
b<-barplot(n_dogs_vax_dist,names=NA,axes=F,space=0.2,col="goldenrod1",border=NA,add=T)
axis(2,at=c(30000,60000,90000),cex.axis=cex.axis,labels=format(c(20000,60000,90000),big.mark=","),padj=1.2)#,at=seq(1,length(years)*12,24),labels=paste(seq(min(years),max(years),2)))
axis(2,cex.axis=cex.axis,at=0,padj=1.2)#,at=seq(1,length(years)*12,24),labels=paste(seq(min(years),max(years),2)))
axis(1,at=b[c(T,F,F,F)],labels=paste(seq(min(years),max(years),4)),cex.axis=cex.axis,padj=-2)
mtext("Dogs",side=2,line=1.6,cex=cex.lab)
mtext("Date",side=1,line=1,cex=cex.lab)
box(bty="l")
legend(-3,90000,legend="B",text.font = 2,bty="n")
legend(1.5,90000,c("Vaccinated this year","Not vaccinated this year"),
       pch=15,col=c("goldenrod1","navy"),bty="n",cex=cex.axis)



## Estimated vax coverage in time
#------------

par(fig=c(0,1,0.63,0.79),new=T)
par(mar=c(1.8,3.5,0.3,3.5))

plot(1:length(vax_dist),vax_dist,type="l",lty=1,col="navy",ylim=c(0,0.8),axes=F,ylab="",xlab="",lwd=2)
axis(1,at=seq(1,length(years)*12,48),labels=paste(seq(min(years),max(years),4)),cex.axis=cex.axis,padj=-2)
axis(2,cex.axis=cex.axis,padj=1.2)
mtext("Rolling vaccination\ncoverage in district",side=2,line=1.6,cex=cex.lab)
mtext("Date",side=1,line=0.9,cex=cex.lab)
box(bty="u")
legend(-20,0.9,legend="C",text.font = 2,bty="n")

par(fig=c(0,1,0.63,0.79),new=T)

plot(1:length(vax_dist),sd_cov,type="l",col="orange",ylim=c(0,0.4),axes=F,ylab="",xlab="",lwd=2,lty=6)
axis(4,cex.axis=cex.axis,padj=-1.5)
mtext("Heterogeneity in\nrolling coverage",side=4,line=2.2,cex=cex.lab)

legend("topright",c("Rolling coverage","Heterogeneity in coverage"),
       lty=c(1,6),col=c("navy","orange"),bty="n",lwd=2,cex=cex.axis)


## Case numbers
#------------

par(fig=c(0,1,0.42,0.64),new=T)
par(mar=c(2.5,2.5,1,0.5))
b<-barplot(c(dog_cases_dist),axes=F,width=1,space=0,col="navy",border="navy",ylim=c(-5,max(dog_cases_dist)))
lines(monthlyBites~b,col="orange",lwd=1)
points(b,rep(-2,length(b)),pch=16,cex=c(0,0.5,0.8)[monthlyDeaths+1],col="firebrick2")
axis(2,at=seq(0,60,20),cex.axis=cex.axis,padj=1)
axis(1,at=seq(1,length(years)*12,24),labels=paste(seq(min(years),max(years),2)),cex.axis=cex.axis,padj=-2)
mtext("Monthly count",side=2,line=1.5,cex=cex.lab)
mtext("Date",side=1,line=1,cex=cex.lab)
box(bty="l")
legend(13.5*12,50,c("Probable dog cases","Probable human exposures"),pch=c(15,NA),pt.cex=c(1.7,NA),col=c("navy","orange"),lty=c(NA,1),lwd=1,bty="n",cex=cex.lab)
legend(13.5*12,70,c("1","2"),title="Human deaths:",pch=16,pt.cex=c(0.4,0.8),col="firebrick2",bty="n",cex=cex.lab,ncol = 2)
legend(-20,75,legend="D",text.font = 2,bty="n")


## Campaign coverage in space
#------------

breaks=seq(0,1,0.01)
colours=colorRampPalette(c("white",brewer.pal(9, "YlOrRd")))(length(breaks)-1)
grid <- raster(extent(SD_vill),crs=SD_vill@proj4string);res(grid) <- 1000;grid[]<-1
nrows<-4
ncols<-6
rows_top <- rep(rev(seq(0,0.42,length.out=nrows+1)[-1]),each=ncols)
rows_bottom <- rep(rev(seq(0,0.42,length.out=nrows+1)[-(nrows+1)]),each=ncols)
cols_start <- rep(seq(0,1,length.out=ncols+1)[-(ncols+1)],nrows)
cols_end <- rep(seq(0,1,length.out=ncols+1)[-1],nrows)
for(i in 1:ncol(campaign_cov_vill)){
  par(fig=c(cols_start[i],cols_end[i],rows_bottom[i],rows_top[i]),mar=c(0,1,0,0),new=T)
  
  plot(SD_vill,col=colours[findInterval(campaign_cov_vill[,i],breaks,all.inside=T)],border=NA)
  plot(SD_outline,add=T)
  points(dog_cases_CT$UTM.Easting[which(dog_cases_CT$year==i+2001)],dog_cases_CT$UTM.Northing[which(dog_cases_CT$year==i+2001)],pch=16,col="blue",cex=cex.pt)

  legend("topright", #121.3,9.8,
         legend=i+2001,text.font=2,
         cex=cex.lab-0.1, bty="n")
  if(i==1){legend(590000,9850000,legend="E",text.font = 2,bty="n",xpd = T)}
}
par(fig=c(cols_start[i+1],cols_end[i+1],rows_bottom[i+1],rows_top[i+1]),mar=c(0,0,0,0),new=T)
plot(grid, 
     breaks=breaks,legend.only=T,col=colours,
     legend.args=list(text="Vaccination\ncampaign\ncoverage", side=4, line=3.3, cex=0.7),
     axis.args=list(at=seq(0,1,0.2) ,labels=seq(0,1,0.2),cex.axis=0.7),
     smallplot=c(0.23,0.28, .07,.69))
par(fig=c(cols_start[i+1],cols_end[i+1],rows_bottom[i+1],rows_top[i+1]),mar=c(0,0,0,0),new=T)
legend("top" ,"Rabid dog",cex=cex.lab-0.1,pt.cex = cex.pt,pch=16,bty="n",col="blue")


dev.off()


## Campaign coverage in district and by village (without cases)
#------------

pdf("Figs/CampaignCoverage.pdf",width=7, height=7)

par(fig=c(0,1,0.7,1),mar=c(2.5,3.5,0.5,1))
plot(percent_dogs_vax/100,type="l",lwd=2,ylim=c(0,1),axes=F,
     ylab="",xlab="")
axis(2,padj=1,cex.axis=cex.axis)
axis(1,at=seq(1,length(vax_dist_annual),2),labels=seq(1,length(vax_dist_annual),2)+2001,padj=-1.5,cex.axis=cex.axis)
mtext("Vaccination campaign\ncoverage in the district",side=2,line=1.5,cex=cex.lab)
mtext("Year",side=1,line=1.2,cex=cex.lab)
legend("topleft","A",text.font=2,cex=1.2,bty="n")
box(bty="l")

breaks=seq(0,1,0.01)
colours=colorRampPalette(c("white",brewer.pal(9, "YlOrRd")))(length(breaks)-1)
grid <- raster(extent(SD_vill),crs=SD_vill@proj4string);res(grid) <- 1000;grid[]<-1
nrows<-4
ncols<-6
rows_top <- rep(rev(seq(0,0.7,length.out=nrows+1)[-1]),each=ncols)
rows_bottom <- rep(rev(seq(0,0.7,length.out=nrows+1)[-(nrows+1)]),each=ncols)
cols_start <- rep(seq(0,1,length.out=ncols+1)[-(ncols+1)],nrows)
cols_end <- rep(seq(0,1,length.out=ncols+1)[-1],nrows)
for(i in 1:ncol(campaign_cov_vill)){
  par(fig=c(cols_start[i],cols_end[i],rows_bottom[i],rows_top[i]),mar=c(0,1,0,0),new=T)
  
  plot(SD_vill,col=colours[findInterval(campaign_cov_vill[,i],breaks,all.inside=T)],border=NA)
  plot(SD_outline,add=T)
  # points(dog_cases_CT$UTM.Easting[which(dog_cases_CT$year==i+2001)],dog_cases_CT$UTM.Northing[which(dog_cases_CT$year==i+2001)],pch=16,col="blue",cex=cex.pt)
  
  if(i==1){legend(597000,9850000,legend="B",text.font = 2,bty="n",xpd = T)}
  legend("topright", #121.3,9.8,
         legend=i+2001,text.font=2,
         cex=cex.lab-0.1, bty="n")
}
par(fig=c(cols_start[i+1],cols_end[i+1],rows_bottom[i+1],rows_top[i+1]),mar=c(0,0,0,0),new=T)
plot(grid, 
     breaks=breaks,legend.only=T,col=colours,
     legend.args=list(text="Vaccination\ncampaign\ncoverage", side=4, line=3.3, cex=0.8),
     axis.args=list(at=seq(0,1,0.2) ,labels=seq(0,1,0.2),cex.axis=0.7,hadj=0.3),
     smallplot=c(0.23,0.28, .15,.85))
par(fig=c(cols_start[i+1],cols_end[i+1],rows_bottom[i+1],rows_top[i+1]),mar=c(0,0,0,0),new=T)

dev.off()



#Mean annual coverage (Fig. S1)
#______________________

pdf("Figs/MeanAnnualCoverage.pdf",width=7, height=8.75)

panel_a_bottom <- 0.75
panel_b_bottom <- 0.5

# time series
set.seed(1)
villages_to_sample <- 8 #nrow(vax_vill_annual)
par(fig=c(0,1,panel_a_bottom,1),mar=c(2.5,3.5,0.5,1))
plot(vax_dist_annual,type="l",col="2",ylim=c(0,1),axes=F,
     ylab="",xlab="")
axis(2,padj=1,cex.axis=cex.axis)
axis(1,at=seq(1,length(vax_dist_annual),2),labels=seq(1,length(vax_dist_annual),2)+2001,padj=-1.5,cex.axis=cex.axis)
mtext("Mean rolling vaccination\ncoverage",side=2,line=1.5,cex=cex.lab)
mtext("Year",side=1,line=1.2,cex=cex.lab)

cols <- brewer.pal(villages_to_sample,"Accent")
samples <-  sample(nrow(vax_vill_annual),villages_to_sample)
for(i in 1:villages_to_sample){lines(vax_vill_annual[samples[i],],col=cols[i],lty=2,lwd=2)}
lines(vax_dist_annual,col=1,lwd=3,ylim=c(0,1))
legend("topright",c("District","Sample villages"),lty=1:2,lwd=3:2,bty="n")
legend("topleft","A",text.font=2,cex=1.2,bty="n")
box(bty="l")

par(fig=c(0,1,panel_b_bottom,panel_a_bottom),mar=c(3,3.5,0.5,1),new=T)
plot(sd_cov_annual,type="l",ylim=c(0,0.2),axes=F,lwd=2,
     ylab="",xlab="")
mtext("Weighted standard deviation\nin rolling coverage",side=2,line=1.5,cex=cex.lab)
mtext("Year",side=1,line=1.2,cex=cex.lab)
axis(2,padj=1,cex.axis=cex.axis)
axis(1,at=seq(1,length(vax_dist_annual),2),labels=seq(1,length(vax_dist_annual),2)+2001,padj=-1.5,cex.axis=cex.axis)
legend("topleft","B",text.font=2,cex=1.2,bty="n")
box(bty="l")

# Maps
cex.pt <- 0.4
breaks=seq(0,1,0.01)
colours=colorRampPalette(c("white",brewer.pal(9, "YlOrRd")))(length(breaks)-1)
grid <- raster(extent(SD_vill),crs=SD_vill@proj4string);res(grid) <- 1000;grid[]<-1
nrows<-4
ncols<-6
rows_top <- rep(rev(seq(0,panel_b_bottom,length.out=nrows+1)[-1]),each=ncols)
rows_bottom <- rep(rev(seq(0,panel_b_bottom,length.out=nrows+1)[-(nrows+1)]),each=ncols)
cols_start <- rep(seq(0,1,length.out=ncols+1)[-(ncols+1)],nrows)
cols_end <- rep(seq(0,1,length.out=ncols+1)[-1],nrows)
for(i in 1:ncol(vax_vill_annual)){
  par(fig=c(cols_start[i],cols_end[i],rows_bottom[i],rows_top[i]),mar=c(0,1,0,0),new=T)
  
  plot(SD_vill,col=colours[findInterval(vax_vill_annual[,i],breaks,all.inside=T)],border=NA)
  plot(SD_outline,add=T)
  # points(dog_cases_CT$UTM.Easting[which(dog_cases_CT$year==i+2001)],dog_cases_CT$UTM.Northing[which(dog_cases_CT$year==i+2001)],pch=16,col="blue3",cex=cex.pt)

  if(i==1){legend(595000,9850000,"C",text.font=2,cex=1.2,bty="n",xpd=T)}
  legend("topright", #121.3,9.8,
         legend=i+2001,text.font=2,
         cex=cex.lab, bty="n")
}
par(fig=c(cols_start[i+1],cols_end[i+1],rows_bottom[i+1],rows_top[i+1]),mar=c(0,0,0,0),new=T)
plot(grid, 
     breaks=breaks,legend.only=T,col=colours,
     legend.args=list(text="Rolling\nvaccination\nCoverage", side=4, line=3.5, cex=0.9),
     axis.args=list(at=seq(0,1,0.2) ,labels=seq(0,1,0.2),cex.axis=0.8,hadj=0.3),
     smallplot=c(0.25,0.29, .15,.85))
par(fig=c(cols_start[i+1],cols_end[i+1],rows_bottom[i+1],rows_top[i+1]),mar=c(0,0,0,0),new=T)
# legend("top" ,"Rabid dog",cex=cex.lab,pt.cex = cex.pt,pch=16,bty="n",col="blue3")


dev.off()


# panel_b_bottom <- 1
# cex.pt <- 0.4
# breaks=seq(0,1,0.01)
# colours=colorRampPalette(c("white",brewer.pal(9, "YlOrRd")))(length(breaks)-1)
# grid <- raster(extent(SD_vill),crs=SD_vill@proj4string);res(grid) <- 1000;grid[]<-1
# nrows<-4
# ncols<-6
# rows_top <- rep(rev(seq(0,panel_b_bottom,length.out=nrows+1)[-1]),each=ncols)
# rows_bottom <- rep(rev(seq(0,panel_b_bottom,length.out=nrows+1)[-(nrows+1)]),each=ncols)
# cols_start <- rep(seq(0,1,length.out=ncols+1)[-(ncols+1)],nrows)
# cols_end <- rep(seq(0,1,length.out=ncols+1)[-1],nrows)
# for(i in 1:ncol(vax_vill_annual)){
#   par(fig=c(cols_start[i],cols_end[i],rows_bottom[i],rows_top[i]),mar=c(0,1,0,0),new=T)
#   values <- ifelse(vax_vill_annual[,i]<0.15,1,0)
#   plot(SD_vill,col=colours[findInterval(values,breaks,all.inside=T)],border=NA)
#   plot(SD_outline,add=T)
#   legend("topright", #121.3,9.8,
#          legend=i+2001,text.font=2,
#          cex=cex.lab, bty="n")
# }



# Cases and vax video
#______________________

breaks=seq(0,1,0.01)
colours=colorRampPalette(c("white",brewer.pal(9, "YlOrRd")))(length(breaks)-1)
grid <- raster(extent(SD_vill),crs=SD_vill@proj4string);res(grid) <- 1000;grid[]<-1
startYear<-min(years)
saveVideo(
  for(i in 1:(length(years)*12)){
    
    par(mar=c(0,0,0,0))
    
    if(i<=(length(years)*12)){
      plot(SD_vill,col=colours[findInterval(vax_vill[,i],breaks,all.inside=T)])
    }else{
      plot(SD_vill,col="grey")
    }
    
    plot(grid, 
         breaks=breaks,legend.only=T, add=T,col=colours,
         legend.args=list(text="Vaccination Coverage", side=4, line=4, cex=2),
         axis.args=list(at=seq(0,1,0.2),cex.axis=1.5),
         smallplot=c(0.80,0.81, .25,.75))
    
    points(dog_cases_CT$UTM.Easting[which(dog_cases_CT$month==i)],dog_cases_CT$UTM.Northing[which(dog_cases_CT$month==i)],pch=16,col="blue",cex=1)
    points(exposures_CT$UTM.Easting[which(exposures_CT$month==i)],exposures_CT$UTM.Northing[which(exposures_CT$month==i)],pch=17,col="magenta3",cex=1)
    points(deaths_CT$UTM.Easting[which(deaths_CT$month==i)],deaths_CT$UTM.Northing[which(deaths_CT$month==i)],pch=17,col="black",cex=2)

    legend("topright" ,c("Rabid dog","Human exposure","Human death"),cex=1.1,pt.cex = c(1,1,2),pch=c(16,17,17),bty="n",col=c("blue","magenta3",1))
    
    legend("topleft", #121.3,9.8,
           legend=paste(month.abb[i - 12*floor((i-1)/12)],floor((i-1)/12)+startYear,sep=" "),
           cex=2, bty="n")
    
    
  },
  video.name = "Figs/Videos/vaccination_cases_exposures_animation_jitter.mp4",
  img.name = "Rplot",
  ffmpeg = ani.options("ffmpeg"),
  interval=0.4,
  ani.width=650,ani.height=480
)


