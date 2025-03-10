#
## Get information about human & dog populations in cells making up Mara grid
#____________________________

rm(list=ls())

library(rgdal)
library(rgeos)
library(lubridate)

options(stringsAsFactors=F) 



## Load data files
#--------------

##Load shapefiles
SD_vill <- readOGR("data/GIS","SD_Villages_2012_From_HHS_250m_Smoothed_UTM") 
SD_vill<-SD_vill[order(SD_vill$Vill_2012),]
SD_vill_2002 <- readOGR("data/GIS","SD_Villages_2002_From_HHS_250m_Smoothed_UTM") 

## Load Serengeti de-identified census data (coordinates and names removed)
Census <- read.csv("data/SDcompiled_deid.csv")

## 2002, 2012 and 2022 NBS census data
wardCensus2002 <- read.csv("data/MaraWardPopNBS2002.csv")
villCensus2012 <- read.csv("data/MaraVillagePopNBS2012.csv")
wardCensus2022 <- read.csv("data/Serengeti_ward_pops_2022.csv")
SDwardCensus2002 <- wardCensus2002[which(wardCensus2002$District=="Serengeti"),]
SDvillCensus2012 <- villCensus2012[which(villCensus2012$District=="Serengeti"),]



## Set dates
#--------------

startYear <- 2000
startYear_analysis <- 2002
endYear <- 2040
years <- length(startYear:endYear)
months <- years*12



## District level human growth rates
#--------------

## Human population size in the district from 2002 NBS
SDpop2002 <- sum(SDwardCensus2002$Population)

## 2012 district population
SDpop2012 <- sum(SDvillCensus2012$Population)

## 2022 district population
SDpop2022 <- sum(wardCensus2022$Population)

## Estimate growth rate 
(Human_growth_2002_2012 <- (SDpop2012/SDpop2002)^(1/(2012-2002))-1)
(Human_growth_2012_2022 <- (SDpop2022/SDpop2012)^(1/(2022-2012))-1)
(Human_growth_2002_2022 <- (SDpop2022/SDpop2002)^(1/(2022-2002))-1)*100
# slightly lower growth rate in 2012-2022 than in 2002-2012



## 2012 Population for villages from NBS
#--------------

## Check village numbers in shapefile against those in census
nrow(SD_vill)
nrow(SDvillCensus2012) #more villages than in shapefile!

## Match shapefile to census dataframe and vice versa
matchVill <- match(SD_vill$Vill_2012,SDvillCensus2012$Village_Match_Groundtruthed_Shapefile)
matchVill2 <- match(SDvillCensus2012$Village_Match_Groundtruthed_Shapefile,SD_vill$Vill_2012)
length(which(is.na(matchVill))) 
SD_vill$Vill_2012[which(is.na(matchVill))]
length(which(is.na(matchVill2))) 
SDvillCensus2012$Village_Match_Groundtruthed_Shapefile[which(is.na(matchVill2))]

## In the groundtruthed shapefile, mugumu and stendi kuu are villages in Mugumu ward. In 2012 census data, the villages are wards

## Serengeti_Nyambureti_Kerukerege occurs in 2012 census but not groundtruthed shapefile. Kerukerege
## seems to have broken off from Maburi. Just put 2012 populations from these two into the Maburi polygon

## Add population data to village shapefile
SD_vill$Humans2012 <- SDvillCensus2012$Population[matchVill]

## Add Kerukerege humans to Maburi 
SD_vill$Humans2012[which(SD_vill$Vill_2012=="Maburi")] <- SD_vill$Humans2012[which(SD_vill$Vill_2012=="Maburi")] + SDvillCensus2012$Population[which(SDvillCensus2012$Village_Match_Groundtruthed_Shapefile=="Kerukerege")]

## Population in Mugumu
SD_vill$Humans2012[which(SD_vill$Vill_2012=="Mugumu")] <- sum(SDvillCensus2012$Population[which(SDvillCensus2012$Ward=="Mugumu")])

## Population in Stendi Kuu
SD_vill$Humans2012[which(SD_vill$Vill_2012=="Stendi Kuu")] <- sum(SDvillCensus2012$Population[which(SDvillCensus2012$Ward=="Stendi Kuu")])



# 2002 government census
#--------------

length(unique(SDwardCensus2002$Ward))
length(unique(SD_vill$Ward_2012))
length(unique(SD_vill_2002$Ward_2002))

# 2002 ward for each 2012 village
coords <- gPointOnSurface(SD_vill, byid = T)
SD_vill$Ward_2002 <- over(coords,SD_vill_2002)$Ward_2002

# Match 2002 ward census data to 2002 village shapefile 
match(SDwardCensus2002$Ward,SD_vill_2002$Ward_2002)
SDwardCensus2002$Ward[which(is.na(match(SDwardCensus2002$Ward,SD_vill_2002$Ward_2002)))]
sort(SD_vill_2002$Ward_2002)
SDwardCensus2002$Ward[which(SDwardCensus2002$Ward=="Kebanchabancha")] <- "Kebanchebanche"
SDwardCensus2002$Ward[which(SDwardCensus2002$Ward=="Mugumu Mjini")] <- "Mugumu Urban"

# Divide 2002 humans in wards among villages based on 2012
SD_vill$propWard02 <- NA
for(ward in unique(SD_vill$Ward_2002)){
  vills <- which(SD_vill$Ward_2002==ward)
  wardPop <- sum(SD_vill$Humans2012[vills])
  SD_vill$propWard02[vills] <- SD_vill$Humans2012[vills]/wardPop
}
SD_vill$wardPop02 <- SDwardCensus2002$Population[match(SD_vill$Ward_2002,SDwardCensus2002$Ward)]
SD_vill$Humans2002 <- SD_vill$wardPop02*SD_vill$propWard02



# 2022 government census
#--------------

# 2012 villages coloured by ward
plot(SD_vill,col=rainbow(length(unique(SD_vill$Ward_2012)))[match(SD_vill$Ward_2012,unique(SD_vill$Ward_2012))])
text(coordinates(SD_vill)[,1],coordinates(SD_vill)[,2],SD_vill$Vill_2012,cex=0.7)

# Assign each village to a 2022 ward
sort(unique(wardCensus2022$Council.Ward))
sort(unique(SD_vill$Ward_2012))
sort(unique(wardCensus2022$Council.Ward))[which(is.na(match(sort(unique(wardCensus2022$Council.Ward)),sort(unique(SD_vill$Ward_2012)))))] # 3 new wards and Kebanchebanche alternative spelling
SD_vill$Ward_2022 <- SD_vill$Ward_2012 # start with 2012 ward names
SD_vill$Ward_2022[which(SD_vill$Vill_2012=="Stendi Kuu")] <- "Stendi Kuu" # Stendi Kuu village now its own ward
SD_vill$Ward_2022[which(SD_vill$Vill_2012=="Matare")] <- "Matare" # new Matare ward
SD_vill$Ward_2022[which(SD_vill$Vill_2012=="Kegonga")] <- "Matare"
SD_vill$Ward_2022[which(SD_vill$Vill_2012=="Kerenero")] <- "Geitasamo" # village moved to neighbouring ward
wardCensus2022$Council.Ward[which(wardCensus2022$Council.Ward=="Kebanchabancha")] <- "Kebanchebanche"
SD_vill$Ward_2022[which(SD_vill$Vill_2012=="Iharara")] <- "Nagusi"
SD_vill$Ward_2022[which(SD_vill$Vill_2012=="Singisi")] <- "Nagusi"
sort(unique(wardCensus2022$Council.Ward))[which(is.na(match(sort(unique(wardCensus2022$Council.Ward)),sort(unique(SD_vill$Ward_2022)))))]
# changes in which villages are assigned to which wards were made using info
# from the postcode directory
# (https://www.tanzaniapostcode.com/search/?keyword=serengeti&region=Mara)

# Geitasamo - Kerenero, Nyameramo, Nyamoko
# Nyamoko - Itununu, Kwitete, Masangura

# Proportion humans in ward in each village
SD_vill$propWard22 <- NA
for(ward in unique(SD_vill$Ward_2022)){
  vills <- which(SD_vill$Ward_2022==ward)
  wardPop <- sum(SD_vill$Humans2012[vills])
  SD_vill$propWard22[vills] <- SD_vill$Humans2012[vills]/wardPop
}
SD_vill$wardPop22 <- wardCensus2022$Population[match(SD_vill$Ward_2022,wardCensus2022$Council.Ward)]
SD_vill$Humans2022 <- SD_vill$wardPop22*SD_vill$propWard22



## Estimate village human Populations through time based on growth rates for each village from NBS census data
#--------------

# Fit growth curve and get population for each village 
(Human_monthly_growth_2002_2022 <- log(SDpop2022/SDpop2002)/((2022-2002)*12)) # initial value for fitting algorithm
humans <- matrix(NA,nrow=nrow(SD_vill),ncol=months)
village_growth_pars <- matrix(NA,nrow=nrow(humans),ncol=2,dimnames=list(NULL,c("a","r")))
for(i in 1:nrow(SD_vill@data)){
  data_i <- data.frame(humans=c(SD_vill$Humans2002[i],SD_vill$Humans2012[i],SD_vill$Humans2022[i]),month=c((2002-2000)*12 + 8,(2012-2000)*12 + 8,(2022-2000)*12 + 8)) # all censuses happened in August of that year
  exp_fit <- nls(humans ~ a*exp(r*month),data=data_i,start=list(a=SD_vill$Humans2002[i],r=Human_monthly_growth_2002_2022))
  village_growth_pars[i,] <- coefficients(exp_fit)
  humans[i,] <- village_growth_pars[i,1]*exp(village_growth_pars[i,2]*c(1:ncol(humans)))
}
par(mar=c(4,4,1,1))
hist(village_growth_pars[,"r"],breaks=20);points(Human_monthly_growth_2002_2022,1,pch=4,col=2,cex=2)
((sum(humans[,12*(length(startYear:2022))+1])/sum(humans[,12*(length(startYear:2002)-1)+1]))^(1/length(2002:2022))-1)*100
((sum(humans[,12*(length(startYear:2022))])/sum(humans[,12*(length(startYear:2002)-1)+1]))^(1/length(2002:2022))-1)*100



## Serengeti census
#--------------

## Proportion dogs that are adults
(prop_adult <- sum(Census$dogs)/sum(Census$dogs + Census$pups))

##HDR for each village
SD_vill$HHS_humans <- c(rowsum(Census$humansTotal, Census$Village))
SD_vill$HHS_dogs <- c(rowsum(Census$dogsTotal, Census$Village))
SD_vill$HDR <- SD_vill$HHS_humans/SD_vill$HHS_dogs

## Overall HDR
sum(Census$humansTotal)/sum(Census$dogsTotal)



## Get dog village populations through time
##--------------

## Get dog population matrices from humans
dogs <- humans/rep(SD_vill$HDR,months)
((sum(dogs[,12*(length(startYear:2022))+1])/sum(dogs[,12*(length(startYear:2002)-1)+1]))^(1/length(2002:2022))-1)*100 # Jan 2002 - Jan 2023



## Save outputs
##--------------

write.table(cbind(SD_vill$Vill_2012,humans),file=paste("output/humanPopulationByVillageMonth_Jan",startYear,"_Dec",endYear,".csv",sep=""), sep=",", row.names = F, col.names = F)
write.table(cbind(SD_vill$Vill_2012,dogs),file=paste("output/dogPopulationByVillageMonth_Jan",startYear,"_Dec",endYear,".csv",sep=""), sep=",", row.names = F, col.names = F)
write.table(cbind(SD_vill$Vill_2012,humans[,((startYear_analysis-startYear)*12+1):months]),file=paste("output/humanPopulationByVillageMonth_Jan",startYear_analysis,"_Dec",endYear,".csv",sep=""), sep=",", row.names = F, col.names = F)
write.table(cbind(SD_vill$Vill_2012,dogs[,((startYear_analysis-startYear)*12+1):months]),file=paste("output/dogPopulationByVillageMonth_Jan",startYear_analysis,"_Dec",endYear,".csv",sep=""), sep=",", row.names = F, col.names = F)
write.table(cbind(SD_vill$Vill_2012,humans[,((startYear_analysis-startYear)*12+1):(12*length(startYear:2022))]),file=paste("output/humanPopulationByVillageMonth_Jan",startYear_analysis,"_Dec2022.csv",sep=""), sep=",", row.names = F, col.names = F)
write.table(cbind(SD_vill$Vill_2012,dogs[,((startYear_analysis-startYear)*12+1):(12*length(startYear:2022))]),file=paste("output/dogPopulationByVillageMonth_Jan",startYear_analysis,"_Dec2022.csv",sep=""), sep=",", row.names = F, col.names = F)
writeOGR(SD_vill, dsn="output/SD_vill", "SD_vill", driver="ESRI Shapefile",overwrite_layer=T)


