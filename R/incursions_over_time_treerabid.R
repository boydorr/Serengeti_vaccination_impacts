rm(list=ls())

library(lubridate)
library(rgdal)
library(viridis)
library(scales)
library(rgeos)



# Prep rabid dog data 
# _________________________

# Import data 
rabid_carnivores <- readRDS(file = "output/clean_bite_data.rda")
n = nrow(rabid_carnivores); n 

# Extract date info
start.date <- as.Date("2002-01-01")
end.date <- as.Date("2022-12-31")
days = 1:length(start.date:end.date); months = month(seq(start.date,end.date,"day")) + (year(seq(start.date,end.date,"day"))-year(start.date))*12
rabid_carnivores$yr = floor((rabid_carnivores$month-1)/12) + year(start.date)

# Add coordinates where missing
no_gps <- which(is.na(rabid_carnivores$UTM.Easting)); length(no_gps); length(no_gps)/n # 43 rabid animals have NO GPS location! (1.2%)
villages_2012 <- readOGR("data/GIS","SD_Villages_2012_From_HHS_250m_Smoothed_UTM") # use village centroid for these
matchVillage <- match(rabid_carnivores$Village[no_gps],villages_2012$Vill_2012)
vill_coords <- coordinates(villages_2012)
rabid_carnivores$UTM.Easting[no_gps] <- vill_coords[matchVillage,1]
rabid_carnivores$UTM.Northing[no_gps] <- vill_coords[matchVillage,2]
plot(villages_2012)
points(rabid_carnivores$UTM.Easting, rabid_carnivores$UTM.Northing, pch=20,col=4) 

# District outline
SD_outline <- gUnaryUnion(villages_2012)# get district outline
SD_outline<-gBuffer(SD_outline,width=1) # get rid of a few tiny holes

# If biter ID does not occur in the ID list, change to 0
rabid_carnivores$Biter.ID[which(!rabid_carnivores$Biter.ID %in% rabid_carnivores$ID & rabid_carnivores$Biter.ID!=0)] <- 0

# Remove duplicates
rabid_carnivores$Duplicated <- duplicated(rabid_carnivores$ID)
dup_ids <- rabid_carnivores[rabid_carnivores$Duplicated,"ID"]
cases = rabid_carnivores[rabid_carnivores$Duplicated==F,] # only keep one record for animals bitten multiple times



# Read in incursions
# _________________________

incursions_99 <- read.csv(file="output/serengeti_incursions_treerabid_prune_99.csv")
incursions_975 <- read.csv(file="output/serengeti_incursions_treerabid_prune_975.csv")
incursions_95 <- read.csv(file="output/serengeti_incursions_treerabid_prune_95.csv")



# Cases and incursions by year
# _________________________

incursions_99$yr <- year(as.Date(incursions_99$date_symptoms))
incursions_975$yr <- year(as.Date(incursions_975$date_symptoms))
incursions_95$yr <- year(as.Date(incursions_95$date_symptoms))

cases_by_year <- hist(cases$yr,breaks=seq(min(cases$yr)-0.5,max(cases$yr)+0.5,1))$counts

incursions_by_year_99 <- hist(incursions_99$yr,breaks=seq(min(cases$yr)-0.5,max(cases$yr)+0.5,1))$counts
incursions_by_year_975 <- hist(incursions_975$yr,breaks=seq(min(cases$yr)-0.5,max(cases$yr)+0.5,1))$counts
incursions_by_year_95 <- hist(incursions_95$yr,breaks=seq(min(cases$yr)-0.5,max(cases$yr)+0.5,1))$counts

prop_incursions_by_year_99 <- incursions_by_year_99/cases_by_year
prop_incursions_by_year_975 <- incursions_by_year_975/cases_by_year
prop_incursions_by_year_95 <- incursions_by_year_95/cases_by_year

# Cut off 2002, because many cases may be wrongly identified as incursions the
# first year
incursions_by_year <- data.frame("year"=2003:2022,
                                 "n_cutoff95"=incursions_by_year_95[-1],
                                 "n_cutoff975"=incursions_by_year_975[-1],
                                 "n_cutoff99"=incursions_by_year_99[-1],
                                 "prop_cutoff95"=round(prop_incursions_by_year_95,2)[-1],
                                 "prop_cutoff975"=round(prop_incursions_by_year_975,2)[-1],
                                 "prop_cutoff99"=round(prop_incursions_by_year_99,2)[-1])

# Summaries for MS
mean(incursions_by_year$n_cutoff95)
mean(incursions_by_year$n_cutoff975)
mean(incursions_by_year$n_cutoff99)
range(incursions_by_year$n_cutoff95)
range(incursions_by_year$n_cutoff975)
range(incursions_by_year$n_cutoff99)
mean(incursions_by_year$prop_cutoff95[1:which(incursions_by_year$year==2017)])
mean(incursions_by_year$prop_cutoff95[which(incursions_by_year$year==2018):nrow(incursions_by_year)])
mean(incursions_by_year$prop_cutoff975[1:which(incursions_by_year$year==2017)])
mean(incursions_by_year$prop_cutoff975[which(incursions_by_year$year==2018):nrow(incursions_by_year)])
mean(incursions_by_year$prop_cutoff99[1:which(incursions_by_year$year==2017)])
mean(incursions_by_year$prop_cutoff99[which(incursions_by_year$year==2018):nrow(incursions_by_year)])
range(incursions_by_year$prop_cutoff95[1:which(incursions_by_year$year==2017)])
range(incursions_by_year$prop_cutoff95[which(incursions_by_year$year==2018):nrow(incursions_by_year)])
range(incursions_by_year$prop_cutoff975[1:which(incursions_by_year$year==2017)])
range(incursions_by_year$prop_cutoff975[which(incursions_by_year$year==2018):nrow(incursions_by_year)])
range(incursions_by_year$prop_cutoff99[1:which(incursions_by_year$year==2017)])
range(incursions_by_year$prop_cutoff99[which(incursions_by_year$year==2018):nrow(incursions_by_year)])
sum(incursions_by_year_99[2:16])/sum(cases_by_year[2:16]) # pre-2018
sum(incursions_by_year_975[2:16])/sum(cases_by_year[2:16]) 
sum(incursions_by_year_95[2:16])/sum(cases_by_year[2:16])
sum(incursions_by_year_99[17:21])/sum(cases_by_year[17:21]) # 2018 onward
sum(incursions_by_year_975[17:21])/sum(cases_by_year[17:21])
sum(incursions_by_year_95[17:21])/sum(cases_by_year[17:21])
sum(incursions_by_year_99[2:17])/sum(cases_by_year[2:17]) # pre-2019
sum(incursions_by_year_975[2:17])/sum(cases_by_year[2:17]) 
sum(incursions_by_year_95[2:17])/sum(cases_by_year[2:17])
sum(incursions_by_year_99[18:21])/sum(cases_by_year[18:21]) # 2019 onward
sum(incursions_by_year_975[18:21])/sum(cases_by_year[18:21])
sum(incursions_by_year_95[18:21])/sum(cases_by_year[18:21])

incursions_by_year$prop_cutoff975
incursions_by_year$prop_cutoff95
incursions_by_year$prop_cutoff99


# Plot incursions and proportion incursions through time
# _________________________

pdf("figs/incursions_treerabid.pdf",7,5)

cex.axis <- 0.7
cex.lab <- 0.8

par(mar=c(2.5,2.5,1,1))
cols <- viridis(3)
cols[3]<- "orange"

par(fig=c(0,0.5,0.5,1))
plot(incursions_by_year_95[-1]~c(2003:2022),type="l",col=cols[1],lwd=1,lty=2,bty="l",ylab="",xlab="",ylim=c(0,30),axes=F)
lines(incursions_by_year_975[-1]~c(2003:2022),col=cols[2],lwd=2,lty=1)
lines(incursions_by_year_99[-1]~c(2003:2022),col=cols[3],lwd=1,lty=2)
axis(2,cex.axis=cex.axis,padj=1)
axis(1,cex.axis=cex.axis,padj=-1.5)
box(bty="l")
mtext("Incursions",side=2,line=1.5,cex=cex.lab)
mtext("Year",side=1,line=1.5,cex=cex.lab)
legend("topleft",legend="A",text.font = 2,bty="n")
legend("topright",title="Pruning threshold:",c("95%","97.5%","99%"),
       title.adj=0,lty=c(2,1,2),title.cex=0.8,cex=0.75,col=cols,lwd=c(1,2,1),bty="n")

par(fig=c(0.5,1,0.5,1),new=T)
plot(prop_incursions_by_year_95[-1]~c(2003:2022),type="l",col=cols[1],lwd=1,lty=2,bty="l",ylab="",xlab="",ylim=c(0,max(prop_incursions_by_year_95)),axes=F)
lines(prop_incursions_by_year_975[-1]~c(2003:2022),col=cols[2],lwd=2,lty=1)
lines(prop_incursions_by_year_99[-1]~c(2003:2022),col=cols[3],lwd=1,lty=2)
axis(2,cex.axis=cex.axis,padj=1)
axis(1,cex.axis=cex.axis,padj=-1)
box(bty="l")
mtext("Proportion incursions",side=2,line=1.5,cex=cex.lab)
mtext("Year",side=1,line=1.5,cex=cex.lab)
legend("topleft",legend="B",text.font = 2,bty="n")
# legend(2010,max(prop_incursions_by_year_95),title="Cut-off:",c("95%","97.5%","95%"),
#        title.adj=0,lty=1:3,title.cex=0.8,cex=0.75,col=cols,lwd=2,bty="n")



# Plot incursions in space
# _________________________


par(mar=c(0,0,2,0))

rabid_carnivores$UTM.Easting_jitter <- jitter(rabid_carnivores$UTM.Easting, amount =1000)
rabid_carnivores$UTM.Northing_jitter <- jitter(rabid_carnivores$UTM.Northing, amount =1000)
incursions_95$x_coord_jitter <- jitter(incursions_95$x_coord, amount =1000)
incursions_95$y_coord_jitter <- jitter(incursions_95$y_coord, amount =1000)
incursions_975$x_coord_jitter <- jitter(incursions_975$x_coord, amount =1000)
incursions_975$y_coord_jitter <- jitter(incursions_975$y_coord, amount =1000)
incursions_99$x_coord_jitter <- jitter(incursions_99$x_coord, amount =1000)
incursions_99$y_coord_jitter <- jitter(incursions_99$y_coord, amount =1000)

par(fig=c(0,0.25,0,0.5),new=T)
plot(villages_2012,border="grey")
plot(SD_outline,add=T)
title("Pruning threshold = 95%",cex.main=0.7,line=0)
points(incursions_95$x_coord_jitter,incursions_95$y_coord_jitter,col=cols[1],pch=16,cex=0.5)
legend(625000,9860000,legend="C",text.font = 2,bty="n",xpd=TRUE)

par(fig=c(0.25,0.5,0,0.5),new=T)
plot(villages_2012,border="grey")
plot(SD_outline,add=T)
title("Pruning threshold = 97.5%",cex.main=0.7,line=0)
points(incursions_975$x_coord_jitter,incursions_975$y_coord_jitter,col=cols[2],pch=16,cex=0.5)

par(fig=c(0.5,0.75,0,0.5),new=T)
plot(villages_2012,border="grey")
plot(SD_outline,add=T)
title("Pruning threshold = 99%",cex.main=0.7,line=0)
points(incursions_99$x_coord_jitter,incursions_99$y_coord_jitter,col=cols[3],pch=16,cex=0.5)

par(fig=c(0.75,1,0,0.5),new=T)
plot(villages_2012,border="grey")
plot(SD_outline,add=T)
title("All cases",cex.main=0.7,line=0)
points(rabid_carnivores$UTM.Easting_jitter,rabid_carnivores$UTM.Northing_jitter,col=alpha("blue3",alpha=1),pch=16,cex=0.5)



dev.off()



