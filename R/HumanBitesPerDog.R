rm(list=ls())

set.seed(0)

endDate <- as.Date("2022-12-31")

library(fitdistrplus)
library(lubridate)


## Read in data
#__________________________

## Human contact tracing
humanCT <- readRDS("output/clean_human_data.rds") 

## Dog cases
rabid_dogs <- read.csv("output/serengeti_rabid_dogs.csv")

# Some checks
humanCTrabid <- humanCT[which(humanCT$Rabid=="Yes" & humanCT$Attacking.species=="Domestic dog"),]
which(!humanCTrabid$Biter.ID%in%rabid_dogs$ID)
table(humanCTrabid$Biter.ID[which(!humanCTrabid$Biter.ID%in%rabid_dogs$ID)],useNA = "ifany") # quite a few humans bitten by rabid dogs not in the animal data
which(humanCT$Biter.ID%in%rabid_dogs$ID & humanCT$Rabid!="Yes") # Not a huge number (14) of cases where the bite is recorded as non-rabid but the dog was suspect
humanCT$Rabid[which(humanCT$Biter.ID%in%rabid_dogs$ID & humanCT$Rabid!="Yes")]



## Get humans bitten by each suspect dog, fit distribution and save
#__________________________

# Humans bitten by each dog
rabid_dogs$humansBitten <- 0
bites_by_ID <- table(humanCT$Biter.ID)
bites_by_ID <- bites_by_ID[-c(1:2)] # remove biter IDs -1 and 0
bites_by_ID <- bites_by_ID[which(as.numeric(names(bites_by_ID))%in%rabid_dogs$ID)] 
rabid_dogs$humansBitten[match(as.numeric(names(bites_by_ID)),rabid_dogs$ID)] <- bites_by_ID
range(rabid_dogs$humansBitten)
mean(rabid_dogs$humansBitten)
length(which(rabid_dogs$humansBitten==0))/nrow(rabid_dogs)
length(which(rabid_dogs$humansBitten==1))/nrow(rabid_dogs)

# Fit negative binomial distribution
human_bites_fit <- fitdist(rabid_dogs$humansBitten,"nbinom")
fitted_values <- dnbinom(0:50,mu=human_bites_fit$estimate["mu"],size=human_bites_fit$estimate["size"])

# Means and CIs
means <- rep(NA,10000)
for(i in 1:10000){
  means[i] <- mean(sample(rabid_dogs$humansBitten,size=nrow(rabid_dogs),replace = T))}
human_bites_CI <- quantile(means,c(0.025,0.975))

# Save outputs
human_bites_per_dog <- rabid_dogs
save(human_bites_per_dog,human_bites_CI,human_bites_fit,file="Output/human_bites_per_dog_outputs.RData")

# Plot 
pdf("Figs/HumanBitesPerDog.pdf",width=4, height=3)
cex.axis <- 0.8
cex.lab <- 0.9
par(mar=c(2.5,2.5,1.5,1))
hist(human_bites_per_dog$humansBitten,breaks=seq(-0.5,max(human_bites_per_dog$humansBitten)+0.5,1),
     ylim=c(0,max(fitted_values,table(human_bites_per_dog$humansBitten)/sum(table(human_bites_per_dog$humansBitten)))),
     border=F,main="",freq=F,xlab="",ylab="",axes=F)
axis(2,cex.axis=cex.axis,padj=1)
axis(1,cex.axis=cex.axis,padj=-1)
mtext("Density",side=2,line=1.5,cex=cex.lab)
mtext("Humans bitten by each suspect dog",side=1,line=1.5,cex=cex.lab)
box(bty="l")
lines(fitted_values~c(0:50),lwd=2,col="navy")
legend("topright",
       c(paste("data (mean=", round(mean(human_bites_per_dog$humansBitten),2),"(95% CI: ", paste0(round(human_bites_CI,2),collapse = "-"),"),\nn=", length(human_bites_per_dog$humansBitten),")",sep=""),
         paste("negative binomial fit\n(mean=", round(human_bites_fit$estimate["mu"],2),", size=",round(human_bites_fit$estimate["size"],2),")",sep="")),
       text.col = c("grey40","navy"),bty="n",y.intersp = 2,cex=cex.lab)
dev.off()
