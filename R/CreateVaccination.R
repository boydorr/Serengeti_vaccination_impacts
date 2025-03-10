rm(list=ls())

library(raster)
library(rgdal)
library(lubridate)
library(RColorBrewer)
library(dplyr)

options(stringsAsFactors=F) 



## Settings
#--------------

## Set dates of interest
startYear <- 2000
startYear_analysis <-2002
endYear <- 2040
years <- length(startYear:endYear)
months <- years*12

## Resolution for grid
res<-1




## Load inputs
#--------------

## Vaccination data (deidentified) from 4 sources
vax <- read.csv(paste("data/VaccinationCleaned.csv",sep=""))
vax_T3 <- read.csv("data/T3vax_serengeti_2020-2022_cleaned.csv")
vax_SHI <- read.csv("data/SHI_vax_aggregated_2021_2022.csv")
vax_thermo <- read.csv("data/thermostability_trial.csv")

## Anna's dog demography data
dog_demog <- read.csv("data/DogdemographydataPlos2016.csv")

## Dog Population data
dogMatVill <- as.matrix(read.csv(paste("output/dogPopulationByVillageMonth_Jan",startYear,"_Dec",endYear,".csv",sep=""),row.names = 1,header=F))

## Load in shapefiles 
vill_2002 <- readOGR("data/GIS","SD_Villages_2002_From_HHS_250m_Lattice_UTM")
vill_2012 <- readOGR("data/GIS","SD_Villages_2012_From_HHS_250m_Lattice_UTM")
vill_2012 <- vill_2012[order(vill_2012$Vill_2012),]
vill_2012_smooth <- readOGR("data/GIS","SD_Villages_2012_From_HHS_250m_Smoothed_UTM")
vill_2012_smooth<-vill_2012_smooth[order(vill_2012_smooth$Vill_2012),]

## Serengeti Census data (deidentified)
Census <- read.csv("data/SDcompiled_deid.csv")

## File for converting 2012 villages to 2002
vill_2002_2012 <- read.csv("data/new2012serengetiVill.csv")



## Prep vaccination data
#___________________________

## Add year and vaccination month to vax
vax$Vaccination.Date <- as.Date(vax$Vaccination.Date)
vax$year <- year(vax$Vaccination.Date)
vax$month <- month(vax$Vaccination.Date)+(year(vax$Vaccination.Date) - (startYear))*12
vax_T3$Vaccination.Date <- as.Date(vax_T3$date.vaccinated)
vax_T3$year <- year(vax_T3$Vaccination.Date)
vax_T3$month <- month(vax_T3$Vaccination.Date)+(year(vax_T3$Vaccination.Date) - (startYear))*12

#Prep SHI data 
vax_SHI <- vax_SHI[,which(!names(vax_SHI)%in%c("cats_vax"))] 
vax_SHI$date <- as.Date(vax_SHI$date)
vax_SHI$year <- year(vax_SHI$date)
vax_SHI$month <- month(vax_SHI$date)+(vax_SHI$year - startYear)*12
vax_SHI$Ward <- vill_2012_smooth$Ward_2012[match(vax_SHI$village,vill_2012_smooth$Vill_2012)]
vax_SHI$District <-"Serengeti"
vax_SHI <- vax_SHI[,c("village","Ward","District","dogs_vax","date","year","month")]
head(vax_SHI)
names(vax_SHI) <- c("Village","Ward","District","Dogs.Vaccinated","Vaccination.Date","year","month")

## Combine T3, SHI, thermostability trial and Wise Monkey data
vax_T3 <- vax_T3[,which(names(vax_T3)!="date.vaccinated")]
names(vax_T3) <- c("Village","Ward","District","Dogs.Vaccinated","Vaccination.Date","year","month")
vax$source <- "Wise Monkey"
vax_T3$source <- "T3"
vax_SHI$source <- "SHI"
vill_in_T3 <- !is.na(match(vill_2012$Vill_2012,vax_T3$Village[which(vax_T3$year==2020)]))
vill_in_WM <- !is.na(match(vill_2012$Vill_2012,vax$Village[which(vax$year==2020)]))
par(mar=c(0,0,0,0))
plot(vill_2012_smooth,col=ifelse(vill_in_T3&vill_in_WM,"purple",ifelse(vill_in_T3,2,ifelse(vill_in_WM,4,"white"))))
legend("topright",c("T3","WM","Both"),fill=c(2,4,"purple"),bty="n")
text(vill_2012,vill_2012$Vill_2012,cex=0.7)
vill_in_T3 <- !is.na(match(vill_2012$Vill_2012,vax_T3$Village[which(vax_T3$year==2021)]))
vill_in_WM <- !is.na(match(vill_2012$Vill_2012,vax$Village[which(vax$year==2021)]))|!is.na(match(vill_2012$Vill_2012,vax_SHI$Village[which(vax_SHI$year==2021)]))
plot(vill_2012_smooth,col=ifelse(vill_in_T3&vill_in_WM,"purple",ifelse(vill_in_T3,2,ifelse(vill_in_WM,4,"white"))))
legend("topright",c("T3","WM","Both"),fill=c(2,4,"purple"),bty="n")
text(vill_2012,vill_2012$Vill_2012,cex=0.7)
vill_in_T3 <- !is.na(match(vill_2012$Vill_2012,vax_T3$Village[which(vax_T3$year==2022)]))
vill_in_WM <- !is.na(match(vill_2012$Vill_2012,vax$Village[which(vax$year==2022)]))|!is.na(match(vill_2012$Vill_2012,vax_SHI$Village[which(vax_SHI$year==2022)]))
plot(vill_2012_smooth,col=ifelse(vill_in_T3&vill_in_WM,"purple",ifelse(vill_in_T3,2,ifelse(vill_in_WM,4,"white"))))
legend("topright",c("T3","WM","Both"),fill=c(2,4,"purple"),bty="n")
text(vill_2012,vill_2012$Vill_2012,cex=0.7)
# plot(vill_2012,col=as.factor(vill_2012$Ward_2012))
vax_thermo$Vaccination.Date <- as.Date(vax_thermo$Vaccination.Date)
vax_thermo$month <- month(vax_thermo$Vaccination.Date)+(vax_thermo$year - startYear)*12
vax <- bind_rows(vax,vax_T3,vax_SHI,vax_thermo)

## Remove vaccinations duplicated between the datasets (facial recognition study
## was carried out in SHI villages, so data is in both Wise Monkey and T3
## datasets)
duplicates<-vax[which(duplicated(vax[,c("Village","Vaccination.Date")]) | duplicated(vax[nrow(vax):1, c("Village","Vaccination.Date")])[nrow(vax):1]),]
duplicates<-duplicates[order(duplicates$Village),]
duplicates_T3 <- duplicates[which(duplicates$source=="T3"),]
vax <- vax[-as.integer(rownames(duplicates_T3)),]

## Want to combine vax data from Kerukerege and Maburi, since they are combined in the shapefile (but not the vaccination data)
vax$Village[which(vax$Village=="Kerukerege")] <- "Maburi"
vax <- aggregate(.~District + Ward + Village + Vaccination.Date + month + year,vax[,c("District","Ward","Village","Dogs.Vaccinated","Vaccination.Date","month","year")],sum,na.rm=T,na.action =NULL)

## Match villages in shapefile to those in vax
vax2012 <- vax[which(vax$year>=2012),]
matchVax <- match(vill_2012$Vill_2012,unique(vax$Village))
matchVax2 <- match(unique(vax$Village),vill_2012$Vill_2012)
vill_2012$Vill_2012[which(is.na(matchVax))]
unique(vax$Village)[which(is.na(matchVax2))]
# Stendi Kuu remains unmatched.  This is not because it wasn't vaccinated, but
# because it was vaccinated simultaneously with Mugumu from the same central
# point

## Split Mugumu vaccinated dogs between Mugumu and Stendi Kuu based on dog population
propDogsMugumu <- dogMatVill["Mugumu",1]/(dogMatVill["Mugumu",1]+dogMatVill["Stendi Kuu",1])
propDogsSK <- 1-propDogsMugumu
for(i in 1:length(which(vax$Village=="Mugumu"))){
  StendiKuu_i <- vax[which(vax$Village=="Mugumu")[i],]
  StendiKuu_i$Village <- "Stendi Kuu"
  StendiKuu_i$Dogs.Vaccinated <- round(propDogsSK*StendiKuu_i$Dogs.Vaccinated)
  vax[which(vax$Village=="Mugumu")[i],"Dogs.Vaccinated"] <- round(propDogsMugumu*vax[which(vax$Village=="Mugumu")[i],"Dogs.Vaccinated"])
  vax <- rbind(vax,StendiKuu_i)
}

## Split Nyamatare vaccinated dogs between Nyamatare and Nyirongo based on dog population from 2020
propDogsNyamatare <- dogMatVill["Nyamatare",1]/(dogMatVill["Nyamatare",1]+dogMatVill["Nyirongo",1])
propDogsNyirongo <- 1-propDogsNyamatare
for(i in 1:length(which(vax$Village=="Nyamatare" & vax$year>=2020))){
  Nyirongo_i <- vax[which(vax$Village=="Nyamatare" & vax$year>=2020)[i],]
  Nyirongo_i$Village <- "Nyirongo"
  Nyirongo_i$Dogs.Vaccinated <- round(propDogsNyirongo*Nyirongo_i$Dogs.Vaccinated)
  vax[which(vax$Village=="Nyamatare" & vax$year>=2020)[i],"Dogs.Vaccinated"] <- round(propDogsNyamatare*vax[which(vax$Village=="Nyamatare" & vax$year>=2020)[i],"Dogs.Vaccinated"])
  vax <- rbind(vax,Nyirongo_i)
}

## There were more villages from 2012 than there were prior to 2012 due to the
## creation of new ones. In most cases new villages were split off from a single
## village.  However, in 2 cases (Tamkeri and Kitarungu) new villages were
## formed from parts split off of two different villages.  Need to split
## vaccinated dogs pre-2012 into the appropriate post-2012 villages in cases
## where the pre-2012 village experienced splitting.

## Sort out the easy cases where a single village broke into 2 or more complete villages
## (and throw out Stendi Kuu because we already fixed that)
easySplits <- vill_2002_2012[which(!vill_2002_2012$village2012 %in% c("Tamkeri","Kitarungu","Sogoti","Stendi Kuu")),]
for(i in 1:length(unique(easySplits$village2002_1))){
  
  ## Find villages that were formed from this village and the pre-2012 vaccination records for the village
  new_villages <- easySplits[which(easySplits$village2002_1==unique(easySplits$village2002_1)[i]),]
  new_villages <- new_villages[order(new_villages$village2012),]
  records <- vax[which(vax$year<2012 & vax$Village==unique(easySplits$village2002_1)[i]),]
  
  ## Get proportions of dogs in each of the 2012 constituent villages
  prop_dogs <- rowsum(Census$dogsTotal[which(Census$Village %in% c(unique(new_villages$village2002_1),new_villages$village2012))],
                        Census$Village[which(Census$Village %in% c(unique(new_villages$village2002_1),new_villages$village2012))])
  prop_dogs <- prop_dogs/sum(prop_dogs)

  ## For each vaccination record, divide up the vaccinated dogs based on prop_dogs
  for(j in 1:nrow(records)){
    newRecords <- records[j,][rep(1,nrow(new_villages)),]
    newRecords$Village <- new_villages$village2012
    newRecords$Dogs.Vaccinated <- round(prop_dogs[which(rownames(prop_dogs)!=unique(new_villages$village2002_1))]*newRecords$Dogs.Vaccinated[1])
    vax[which(vax$year<2012 & vax$Village==unique(easySplits$village2002_1)[i])[j],"Dogs.Vaccinated"] <- 
      round(prop_dogs[which(rownames(prop_dogs)==unique(new_villages$village2002_1))]*records[j,"Dogs.Vaccinated"])
    vax <- rbind(vax,newRecords)
  }
}


## More complex cases (Tamkeri, Kitarungu, Sogoti, and their 4 precursers)
records_Kebanchebanche <- vax[which(vax$year<2012 & vax$Village=="Kebanchebanche"),]
records_Mbalibali <- vax[which(vax$year<2012 & vax$Village=="Mbalibali"),]
records_Nyansurura <- vax[which(vax$year<2012 & vax$Village=="Nyansurura"),]
records_Nyamburi <- vax[which(vax$year<2012 & vax$Village=="Nyamburi"),]
props_Kebanchebanche <- rowsum(Census$dogsTotal[which(Census$Village.2002 =="Kebanchebanche")],
                               Census$Village[which(Census$Village.2002 =="Kebanchebanche")])
props_Kebanchebanche <- props_Kebanchebanche/sum(props_Kebanchebanche)
props_Mbalibali <- rowsum(Census$dogsTotal[which(Census$Village.2002 =="Mbalibali")],
                               Census$Village[which(Census$Village.2002 =="Mbalibali")])
props_Mbalibali <- props_Mbalibali/sum(props_Mbalibali)
props_Nyansurura <- rowsum(Census$dogsTotal[which(Census$Village.2002 =="Nyansurura")],
                               Census$Village[which(Census$Village.2002 =="Nyansurura")])
props_Nyansurura <- props_Nyansurura/sum(props_Nyansurura)
props_Nyamburi <- rowsum(Census$dogsTotal[which(Census$Village.2002 =="Nyamburi")],
                               Census$Village[which(Census$Village.2002 =="Nyamburi")])
props_Nyamburi <- props_Nyamburi/sum(props_Nyamburi)


## Sogoti
for(i in 1:nrow(records_Kebanchebanche)){
  Sogoti_i <- records_Kebanchebanche[i,]
    Sogoti_i$Village <- "Sogoti"
    Sogoti_i$Dogs.Vaccinated <- round(props_Kebanchebanche[which(rownames(props_Kebanchebanche)=="Sogoti")]*
                                        Sogoti_i$Dogs.Vaccinated)
    vax <- rbind(vax,Sogoti_i)}


## Kitarungu
for(i in 1:nrow(records_Kebanchebanche)){
  Kitarungu_i <- records_Kebanchebanche[i,]
  Kitarungu_i$Village <- "Kitarungu"
  Kitarungu_i$Dogs.Vaccinated <- round(props_Kebanchebanche[which(rownames(props_Kebanchebanche)=="Kitarungu")]*
                                         Kitarungu_i$Dogs.Vaccinated)
  vax <- rbind(vax,Kitarungu_i)}
for(i in 1:nrow(records_Nyansurura)){
  Kitarungu_i <- records_Nyansurura[i,]
  Kitarungu_i$Village <- "Kitarungu"
  Kitarungu_i$Dogs.Vaccinated <- round(props_Nyansurura[which(rownames(props_Nyansurura)=="Kitarungu")]*
                                         Kitarungu_i$Dogs.Vaccinated)
  vax <- rbind(vax,Kitarungu_i)}


## Tamkeri
for(i in 1:nrow(records_Nyamburi)){
  Tamkeri_i <- records_Nyamburi[i,]
  Tamkeri_i$Village <- "Tamkeri"
  Tamkeri_i$Dogs.Vaccinated <- round(props_Nyamburi[which(rownames(props_Nyamburi)=="Tamkeri")]*
                                         Tamkeri_i$Dogs.Vaccinated)
  vax <- rbind(vax,Tamkeri_i)}
for(i in 1:nrow(records_Mbalibali)){
  Tamkeri_i <- records_Mbalibali[i,]
  Tamkeri_i$Village <- "Tamkeri"
  Tamkeri_i$Dogs.Vaccinated <- round(props_Mbalibali[which(rownames(props_Mbalibali)=="Tamkeri")]*
                                         Tamkeri_i$Dogs.Vaccinated)
  vax <- rbind(vax,Tamkeri_i)}


## Kebanchebanche
vax$Dogs.Vaccinated[which(vax$year<2012 & vax$Village=="Kebanchebanche")] <-
  round(vax$Dogs.Vaccinated[which(vax$year<2012 & vax$Village=="Kebanchebanche")]*
          props_Kebanchebanche[which(rownames(props_Kebanchebanche)=="Kebanchebanche")])

## Mbalibali
vax$Dogs.Vaccinated[which(vax$year<2012 & vax$Village=="Mbalibali")] <-
  round(vax$Dogs.Vaccinated[which(vax$year<2012 & vax$Village=="Mbalibali")]*
          props_Mbalibali[which(rownames(props_Mbalibali)=="Mbalibali")])

## Nyansurura
vax$Dogs.Vaccinated[which(vax$year<2012 & vax$Village=="Nyansurura")] <-
  round(vax$Dogs.Vaccinated[which(vax$year<2012 & vax$Village=="Nyansurura")]*
          props_Nyansurura[which(rownames(props_Nyansurura)=="Nyansurura")])

## Nyamburi
vax$Dogs.Vaccinated[which(vax$year<2012 & vax$Village=="Nyamburi")] <-
  round(vax$Dogs.Vaccinated[which(vax$year<2012 & vax$Village=="Nyamburi")]*
          props_Nyamburi[which(rownames(props_Nyamburi)=="Nyamburi")])

## Aggregate records from same village in same month
vax_agg <- aggregate(.~month + Village + year, vax[,c("Village","month","year","Dogs.Vaccinated")],sum)





## Create matrix of vaccinated dogs per village per month and save
#___________________________

## Create matrix to hold vaccination events  by village
vaxVill <- matrix(0, ncol=months, nrow=nrow(vill_2012))

## Add dogs vaccinated in each month to vaxVill
vaxVill[cbind(match(vax_agg$Village,unique(vill_2012$Vill_2012)),vax_agg$month)] <- vax_agg$Dogs.Vaccinated
vaxVill_no2001campaign <- vaxVill

## Add 2001 vaccination
## 73.7% of dogs were vaccinated between May 2000 and February 2001
wards_vaxed_2001 <- c("Issenye","Natta","Ikoma","Sedeco","Kyambahi","Manchira","Rigicha","Kisangura","Nyambureti","Mosongo","Uwanja wa Ndege","Morotonga","Mugumu","Mbalibali","Machochwe","Geitasamo","Rung'abure","Nyansurura","Nyamoko","Kebanchebanche") #wards I think were part of SD at this time based on map from Cleaveland 2003 
plot(vill_2012_smooth)
plot(vill_2012_smooth[which(vill_2012_smooth$Ward_2012%in%wards_vaxed_2001),],add=T,col=3)
vaxVill[which(vill_2012_smooth$Ward_2012%in%wards_vaxed_2001),5:14] <- round((0.737*rowMeans(dogMatVill[which(vill_2012_smooth$Ward_2012%in%wards_vaxed_2001),5:14]))/length(5:14))

## Save dogs vaccinated matrix
dogMatVill <- dogMatVill[order(match(rownames(dogMatVill),vill_2012$Vill_2012)),]
write.table(cbind(vill_2012$Vill_2012,vaxVill),paste("output/dogsVaccinatedByVillageMonth_Jan",startYear,"_Dec",endYear,".csv",sep=""),
            col.names = F,row.names = F,sep=",")

## Vaccination by year
vaxVillYear <- matrix(0,nrow=nrow(vill_2012),ncol=ncol(vaxVill)/12)
for(i in 1:ncol(vaxVillYear)){
  vaxVillYear[,i]<-rowSums(vaxVill[,(1:12)+12*(i-1)])
}
write.table(vaxVillYear[,3:23],paste("output/dogsVaccinatedByVillageYear.csv",sep=""),
            col.names = 2002:2022,row.names = F,sep=",")
vaxVillYearBern <- vaxVillYear[,3:23]
vaxVillYearBern[which(vaxVillYearBern>0)]<-1
par(mfrow=c(3,3))
for(i in 1:ncol(vaxVillYearBern)){
  plot(vill_2012,main=2002+i-1,col=vaxVillYearBern[,i])
}
write.table(vaxVillYearBern,paste("output/VillagesVaccinatedByYear.csv",sep=""),
            col.names = 2002:2022,row.names = F,sep=",")




## Campaign coverage by month
#___________________________

## Check for coverages of >1 (& for NAs)
vc <- vaxVill/dogMatVill
which(is.na(vc))
length(which(vc>1))/length(which(vc>0))
cbind(which(vc>1,arr.ind = T),vc[which(vc>1,arr.ind = T)]) # Matare accounts for 50% of coverages >1
# Some cases but not excessive 

## Limit coverage 
a<-6 # 100% gets corrected to 89%, 70% to 69%, but 60% and under is mostly unchanged
vc_adjust <- vc/(1+vc^a)^(1/a)
max(vc_adjust)

## Save monthly matrix with just campaign coverage
write.table(cbind(vill_2012$Vill_2012,vc_adjust[,((startYear_analysis-startYear)*12+1):ncol(vc)]),paste("output/vaccinationCampaignCoverageByVillageMonth_Jan",startYear_analysis,"_Dec",endYear,".csv",sep=""),
            col.names = F,row.names = F,sep=",")


# Explore what happened with the missed SHI campaign in 2021
SHI_villages <- c(unique(vax_SHI$Village[which(year(vax_SHI$Vaccination.Date)<2023)]),"Stendi Kuu","Kitembere")
missed_2021 <- vill_2012_smooth$Vill_2012[which(vaxVillYearBern[,2021-2002+1]==0)]
missed_2021_SHI <- SHI_villages[which(SHI_villages%in%missed_2021)]
length(SHI_villages)/nrow(vill_2012_smooth)
missed_2021 <- which(vaxVillYearBern[,2021-2002+1]==0)
gap <- matrix(NA,ncol=4,nrow=length(missed_2021))
for(i in 1: length(missed_2021)){
  non_zero <- which(vaxVill[missed_2021[i],]>0)
  gap[i,1] <- vill_2012_smooth$Vill_2012[missed_2021[i]]
  gap[i,2] <- max(non_zero[which(non_zero<((2021-2000)*12+1))])
  gap[i,3] <- min(non_zero[which(non_zero>((2021-2000+1)*12))])
  gap[i,4] <- as.numeric(gap[i,3]) - as.numeric(gap[i,2])
}
gap 
table(as.numeric(gap[,4])) # only village where the gap was >18 month was Nyamatoke (which was not added to SHI program until 2023)
gap_shi <- gap[which(gap[,1]%in%SHI_villages),]
table(as.numeric(gap_shi[,4]))

# Explore what happened with the missed northwest campaign in 2018
missed_2018 <- which(vaxVillYearBern[,2018-2002+1]==0)
gap <- matrix(NA,ncol=4,nrow=length(missed_2018))
for(i in 1: length(missed_2018)){
  non_zero <- which(vaxVill[missed_2018[i],]>0)
  gap[i,1] <- vill_2012_smooth$Vill_2012[missed_2018[i]]
  gap[i,2] <- max(non_zero[which(non_zero<((2018-2000)*12+1))])
  gap[i,3] <- min(non_zero[which(non_zero>((2018-2000+1)*12))])
  gap[i,4] <- as.numeric(gap[i,3]) - as.numeric(gap[i,2])
}
gap 
table(as.numeric(gap[,4])) # all >1.5 years




## Campaign coverage by year (assumes no revax of same dogs within a year)
#___________________________

## Sum coverages by year 
vcVillYear <- matrix(0,nrow=nrow(vc),ncol=ncol(vc)/12)
for(i in 1:ncol(vcVillYear)){
  vcVillYear[,i]<-rowSums(vc[,(1:12)+12*(i-1)])
}
vcVillYear <- vcVillYear[,(2002-startYear+1):(length(startYear:2022))]
rownames(vcVillYear) <- vill_2012$Vill_2012
colnames(vcVillYear) <- 2002:2022
hist(vcVillYear)
max(vcVillYear)
length(which(vcVillYear>1))
length(which(vcVillYear>1))/length(which(vcVillYear>0))
1-length(which(vcVillYear>0.7))/length(which(vcVillYear>0))
tail(sort(vcVillYear))
which(vcVillYear>1,arr.ind = T)
vcVillYear <- vcVillYear/(1+vcVillYear^a)^(1/a)
max(vcVillYear)
write.table(vcVillYear,"output/vcVillByYear.csv",col.names = F,row.names = F,sep=",")

## Get annual vaccination coverage by district
vcDist <- colSums(vaxVill)/colSums(dogMatVill) 
vcDistYear <- rep(0,ncol(vc)/12)
for(i in 1:length(vcDistYear)){
  vcDistYear[i]<-sum(vcDist[(1:12)+12*(i-1)])
}
max(vcDistYear)
write.table(vcDistYear[(2002-startYear+1):length(startYear:2022)],"output/vcDistByYear.csv",row.names = F,col.names = F,sep=",")



## Estimate dog death rate
#___________________________

## dog demography data was collected annually for each participating household during 2010-2013 (though not all
## households/dogs will have data for all four years)
dog_demog$deathdate <- as.Date(dog_demog$deathdate,format="%d/%m/%Y")
dog_demog$birthdate <- as.Date(gsub(" ", "", dog_demog$birthdate),format="%m/%d/%Y")
dog_demog$encdate <- as.Date(gsub(" ", "", dog_demog$encdate, fixed = TRUE),format="%m/%d/%Y")
dog_demog$year_recorded <- year(dog_demog$encdate)
which(dog_demog$year!=dog_demog$year_recorded)

# Proportion of dogs that are pups in each year
dog_demog %>% group_by(year_recorded) %>% summarise(length(which(age<3))/n())
dog_demog %>% group_by(year_recorded) %>% summarise(length(which(age<=3))/n())
sum(Census$pups)/sum(Census$dogsTotal) 
# Looks ok in first couple years, if using 3 month olds but otherwise not

##Throw out dogs that were only observed in one year
single_obs <- names(table(dog_demog$dogid))[which(table(dog_demog$dogid)==1)]
length(which(dog_demog$age[which(dog_demog$dogid%in%single_obs)]<3))
length(which(dog_demog$age[which(dog_demog$dogid%in%single_obs)]>=3))
length(which(dog_demog$age[which(dog_demog$dogid%in%single_obs)]<3))/length(single_obs)#Proportion of single observations that are pups is higher than the proportion of dogs that are pups
table(dog_demog$year[which(dog_demog$dogid%in%single_obs)]) # all but one in the last year (so don't seem to be losing dogs that changed household)
dog_demog <- dog_demog[-which(is.element(dog_demog$dogid, single_obs)),] 


# cases where the dog died but there was no visit in the prior year
dog_demog <- 
  dog_demog %>% group_by(dogid) %>% 
  mutate(death_year=min(year_recorded[which(dead==1)]),death_year=ifelse(!is.infinite(death_year),death_year,NA),
         visit_before_death=(death_year-1)%in%year_recorded) %>% 
  ungroup() 
length(which(dog_demog$year_recorded==dog_demog$death_year & dog_demog$visit_before_death==FALSE)) #only one case!
dog_demog$deathdate[which(dog_demog$year_recorded==dog_demog$death_year & dog_demog$visit_before_death==FALSE)] # and it does not have a death date
dog_demog$dogid[which(dog_demog$year_recorded==dog_demog$death_year & dog_demog$visit_before_death==FALSE)] 
dog_demog$encdate[which(dog_demog$dogid=="ums082")] 
# can just ignore - the loss of possibly one survival year will make basically
# no difference to the death rate estimated below

## Get subsets of pup and adult survivals
unique(dog_demog%>%group_by(dogid)%>%mutate(always_dead=any(dead==1)&!any(alive==1))%>%pull(always_dead))
length(which(is.na(dog_demog$age))) # lots - primarily because dead dogs are typically given NA as age
table(dog_demog$age[which(dog_demog$dead==1)],useNA = "ifany")
unknown_age <- dog_demog[which(is.na(dog_demog$age)& dog_demog$dead==0),]
(unknown_age$encdate-unknown_age$birthdate)/30 # all >3months
dog_demog_adult <- dog_demog[which(dog_demog$age>=3|is.na(dog_demog$age)),]
single_obs <- names(table(dog_demog_adult$dogid))[which(table(dog_demog_adult$dogid)==1)]
dog_demog_adult <- dog_demog_adult[-which(is.element(dog_demog_adult$dogid, single_obs)),] 
dog_demog_adult <- dog_demog_adult %>% group_by(dogid) %>% mutate(always_dead=length(unique(dead))==1 & unique(dead)[1]==1) %>% filter(always_dead==FALSE)
dog_demog_pup <- dog_demog %>% group_by(dogid) %>% mutate(pup=any(age<3)) %>% filter(pup) %>% slice_min(order_by=year,n=2) %>% filter(diff(year)==1)
unique(dog_demog %>% group_by(dogid) %>% mutate(households=length(unique(hhid))) %>% pull(households)) # no household movement (no cases where pup moved hh? Could be missing survivals if these are missed - check with Anna)

## number of cases where a dog died during the study period
deaths <- length(unique(dog_demog$dogid[which(dog_demog$dead==1)])) # 1559 deaths of study dogs
deaths_adults <- length(unique(dog_demog_adult$dogid[which(dog_demog_adult$dead==1)])) #992
deaths_pups <- length(unique(dog_demog_pup$dogid[which(dog_demog_pup$dead==1)])) #567

## number of times dogs survived between years
survivals <- survivals_adults <- survivals_pups <- 0
for(i in 1:length(unique(dog_demog$dogid))){
  obs <- which(dog_demog$dogid==unique(dog_demog$dogid)[i])
  survivals <- survivals + sum(diff(dog_demog$year[obs])[which(dog_demog$alive[obs[-1]]==1)]) 
}
for(i in 1:length(unique(dog_demog_adult$dogid))){
  obs <- which(dog_demog_adult$dogid==unique(dog_demog_adult$dogid)[i])
  survivals_adults <- survivals_adults + sum(diff(dog_demog_adult$year[obs])[which(dog_demog_adult$alive[obs[-1]]==1)]) 
}
for(i in 1:length(unique(dog_demog_pup$dogid))){
  obs <- which(dog_demog_pup$dogid==unique(dog_demog_pup$dogid)[i])
  survivals_pups <- survivals_pups + sum(diff(dog_demog_pup$year[obs])[which(dog_demog_pup$alive[obs[-1]]==1)]) 
  cat(diff(dog_demog_pup$year[obs]))
}
survivals # 1918 cases where dog survived the year
survivals_adults # 1834 cases where adult dog survived the year
survivals_pups # 84 cases where pup survived the year

## estimate death rate
(annual_death_rate <- deaths/(survivals+deaths)) #proportion of dogs that die in a given year
(annual_death_rate_adult <- deaths_adults/(survivals_adults+deaths_adults)) # 0.35
(annual_death_rate_pup <- deaths_pups/(survivals_pups+deaths_pups))  # pup deaths very high 

(survivals_pups+deaths_pups)/(survivals_pups+deaths_pups+survivals_adults+deaths_adults)
annual_death_rate_pup*0.17+annual_death_rate_adult*(1-0.17)



## Estimate monthly coverage/immunity 
#___________________________

## Village-level 
#----------

## Duration of vaccine-induced rabies immunity
vaccine_duration <- 3 # years

## Monthly rate of vaccination waning as a consequence of replacement of deaths and vaccine duration
# lambda <- exp(-( - log(1-annual_death_rate))/12)
lambda <- exp(-(1/vaccine_duration - log(1-annual_death_rate))/12)
lambda_turnover <- exp(-( - log(1-annual_death_rate))/12)

## Create matrix to hold vaccination coverages by village 
vc_waning <- vaxVill
immune_waning <- vaxVill
years <- rep(1:(ncol(vaxVill)/12),each=12) # group columns by year
years[1:12] <- 2# make 2000 and 2001 the same year
previous_years_vax <- rep(0,nrow(vc_waning)) # initialise vector describing number of vaccinated dogs that were not vaccinated in the current year
previous_years_immune <- rep(0,nrow(vc_waning)) # initialise vector describing number of immune dogs that were not vaccinated in the current year
for(i in 2:ncol(vc_waning)){ # each month
  if(years[i]>years[i-1]){ # Jan
    previous_years_immune <- lambda*immune_waning[,i-1] # all the currently immune dogs were vaxed in previous years (with some waning from last month given by lambda)
    previous_years_vax <- lambda_turnover*vc_waning[,i-1] # all the currently vaxed dogs were vaxed in previous years (with some waning from last month given by lambda)
    vax_this_month <- vc_waning[,i] # number of dogs to be vaccinated this month
    immune_waning[,i] <- previous_years_immune + vax_this_month*(1-(previous_years_immune/dogMatVill[,i]))
    immune_waning[,i] <- ((immune_waning[,i]/dogMatVill[,i])/(1+(immune_waning[,i]/dogMatVill[,i])^a)^(1/a))*dogMatVill[,i]
    vc_waning[,i] <- previous_years_vax + vax_this_month*(1-(previous_years_vax/dogMatVill[,i])) # randomly vaccinate amongst dogs 
    vc_waning[,i] <- ((vc_waning[,i]/dogMatVill[,i])/(1+(vc_waning[,i]/dogMatVill[,i])^a)^(1/a))*dogMatVill[,i]
    # number of vaccinated dogs is the same as last month +waning +any new
    # vaccinations applied to susceptibles.  Can't exceed dog population.
    previous_years_vax <- pmax(0, previous_years_vax - vax_this_month*(previous_years_vax/dogMatVill[,i]))  
    previous_years_immune <- pmax(0,previous_years_immune-vax_this_month*(previous_years_immune/dogMatVill[,i]))
    
  }else{ # Feb-Dec
    previous_years_vax <- lambda_turnover*previous_years_vax
    previous_years_immune <- lambda*previous_years_immune 
    vax_this_month <- vc_waning[,i] # number of dogs to be vaccinated this month
    immune_waning[,i] <- lambda*immune_waning[,i-1] + vax_this_month*(1-(previous_years_immune/(dogMatVill[,i]-(lambda_turnover*vc_waning[,i-1]-previous_years_vax))))
    immune_waning[,i] <- ((immune_waning[,i]/dogMatVill[,i])/(1+(immune_waning[,i]/dogMatVill[,i])^a)^(1/a))*dogMatVill[,i]
    vc_waning[,i] <- lambda_turnover*vc_waning[,i-1] + vax_this_month*(1-(previous_years_vax/(dogMatVill[,i]-(lambda_turnover*vc_waning[,i-1]-previous_years_vax)))) # randomly vaccinate amongst dogs NOT ALREADY VAXED THIS YEAR
    vc_waning[,i] <- ((vc_waning[,i]/dogMatVill[,i])/(1+(vc_waning[,i]/dogMatVill[,i])^a)^(1/a))*dogMatVill[,i]
    previous_years_vax <- pmax(0, previous_years_vax - vax_this_month*(previous_years_vax/(dogMatVill[,i]-(lambda_turnover*vc_waning[,i-1]-previous_years_vax)))) 
    previous_years_immune <- pmax(0,previous_years_immune-vax_this_month*(previous_years_immune/(dogMatVill[,i]-(lambda_turnover*vc_waning[,i-1]-previous_years_vax))))
  }
}

vc_waning<-vc_waning/dogMatVill
which((vc_waning)>1)
which(is.na(vc_waning))
mean(vc_waning[,25]) 

immune_waning<-immune_waning/dogMatVill
which((immune_waning)>1)
which(is.na(immune_waning))
mean(immune_waning[,25]) 

## Save vaccination coverage matrix
write.table(cbind(vill_2012$Vill_2012,vc_waning[,((startYear_analysis-startYear)*12+1):ncol(vc_waning)]),paste("output/vaccinationCoverageByVillageMonth_Jan",startYear_analysis,"_Dec",endYear,".csv",sep=""),
            col.names = F,row.names = F,sep=",")
write.table(cbind(vill_2012$Vill_2012,vc_waning),paste("output/vaccinationCoverageByVillageMonth_Jan",startYear,"_Dec",endYear,".csv",sep=""),
            col.names = F,row.names = F,sep=",")
write.table(cbind(vill_2012$Vill_2012,vc_waning[,((startYear_analysis-startYear)*12+1):(12*length(startYear:2022))]),paste("output/vaccinationCoverageByVillageMonth_Jan",startYear_analysis,"_Dec2022.csv",sep=""),
            col.names = F,row.names = F,sep=",")

## Save immunity matrix
write.table(cbind(vill_2012$Vill_2012,immune_waning[,((startYear_analysis-startYear)*12+1):ncol(immune_waning)]),paste("output/immunityByVillageMonth_Jan",startYear_analysis,"_Dec",endYear,".csv",sep=""),
            col.names = F,row.names = F,sep=",")
write.table(cbind(vill_2012$Vill_2012,immune_waning),paste("output/immunityByVillageMonth_Jan",startYear,"_Dec",endYear,".csv",sep=""),
            col.names = F,row.names = F,sep=",")
write.table(cbind(vill_2012$Vill_2012,immune_waning[,((startYear_analysis-startYear)*12+1):(12*length(startYear:2022))]),paste("output/immunityByVillageMonth_Jan",startYear_analysis,"_Dec2022.csv",sep=""),
            col.names = F,row.names = F,sep=",")



## District-level
#----------

par(mar=c(5,5,1,1))
par(mfrow=c(1,1))
district_ts <- data.frame("dogsVax"=colSums(vc_waning*dogMatVill),
                          "dogsImmune"=colSums(immune_waning*dogMatVill))
district_ts$vc <- district_ts$dogsVax/colSums(dogMatVill)
district_ts$immunity <- district_ts$dogsImmune/colSums(dogMatVill)
plot(district_ts$vc,type="l")
lines(district_ts$immunity,col="red")

## Save district-level info
write.table(district_ts,paste("output/districtVaccinationCoverage_Jan",startYear,"_Dec",endYear,".csv",sep=""),
            col.names = T,row.names = F,sep=",")
write.table(district_ts[((startYear_analysis-startYear)*12+1):nrow(district_ts),],paste("output/districtVaccinationCoverage_Jan",startYear_analysis,"_Dec",endYear,".csv",sep=""),
            col.names = T,row.names = F,sep=",")
write.table(district_ts[((startYear_analysis-startYear)*12+1):(12*length(startYear:2022)),c(3,4)],paste("output/districtVaccinationCoverage_Jan",startYear_analysis,"_Dec2022.csv",sep=""),
            col.names = T,row.names = F,sep=",")


