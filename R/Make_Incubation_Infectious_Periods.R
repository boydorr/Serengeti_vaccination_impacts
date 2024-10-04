rm(list=ls())
library(ggplot2)
library(fitdistrplus)
require(MASS)
library(DHARMa)

source("R/Functions/computePeriod.R")

set.seed(0)

# This file is to be run after Data_Cleaning.R
# To do this, we need to use the augmented rabid_carnivores dataset
rabid_carnivores <- readRDS(file = "output/clean_bite_data.rda"); nrow(rabid_carnivores) 

# Identify records of animals bitten by multiple other rabid animals (and so we don't know true date bitten)
MultiBite <- which(rabid_carnivores$ID %in% rabid_carnivores$ID[which(duplicated(rabid_carnivores$ID))])



## Incubation Period
#_______________________

# Convert incubation periods into units of days
# See explanation of algorithm in computePeriod.
rabid_carnivores$Incubation.Period.Clean <- computePeriod(rabid_carnivores$Incubation.period,
                              rabid_carnivores$Incubation.period.units,
                              rabid_carnivores$Date.bitten,
                              rabid_carnivores$Symptoms.started,
                              rabid_carnivores$Date.bitten.uncertainty,
                              rabid_carnivores$Symptoms.started.accuracy)
rabid_carnivores$Incubation.Period.Clean[MultiBite] <- NA # discard incubation periods for dogs that have multiple possible values
summary(rabid_carnivores$Incubation.Period.Clean)

## Two cases with incubation period recorded as zero - change to 3 
which(rabid_carnivores$Incubation.Period.Clean==0) 
rabid_carnivores$Incubation.Period.Clean[which(rabid_carnivores$Incubation.Period.Clean==0)]<-3 # change 0 period to 3


## Fit incubation period distributions
inc_idx <- which(rabid_carnivores$Incubation.Period.Clean>=1); length(inc_idx) # 1212
(inc_gamma = fitdist(rabid_carnivores$Incubation.Period.Clean[inc_idx], "gamma"))
(inc_lnorm = fitdist(rabid_carnivores$Incubation.Period.Clean[inc_idx], "lnorm"))
(inc_weibull = fitdist(rabid_carnivores$Incubation.Period.Clean[inc_idx], "weibull"))
round(c(inc_gamma$aic, inc_lnorm$aic, inc_weibull$aic)) # lognormal is best
boot_inc_gamma <- summary(bootdist(inc_gamma, niter = 1001))
boot_inc_lnorm <- summary(bootdist(inc_lnorm, niter = 1001))
boot_inc_weibull <- summary(bootdist(inc_weibull, niter = 1001))
inc_MS_summary <- rbind(c(paste0(round(inc_gamma$estimate,2)," (",round(boot_inc_gamma$CI[,"2.5%"],2),"-",round(boot_inc_gamma$CI[,"97.5%"],2),")"),round(inc_gamma$aic),round(qgamma(c(0.95,0.975,0.99),shape=inc_gamma$estimate["shape"],rate=inc_gamma$estimate["rate"]),1)),
                        c(paste0(round(inc_lnorm$estimate,2)," (",round(boot_inc_lnorm$CI[,"2.5%"],2),"-",round(boot_inc_lnorm$CI[,"97.5%"],2),")"),round(inc_lnorm$aic),round(qlnorm(c(0.95,0.975,0.99),meanlog=inc_lnorm$estimate["meanlog"],sdlog=inc_lnorm$estimate["sdlog"]),1)),
                        c(paste0(round(inc_weibull$estimate,2)," (",round(boot_inc_weibull$CI[,"2.5%"],2),"-",round(boot_inc_weibull$CI[,"97.5%"],2),")"),round(inc_weibull$aic),round(qweibull(c(0.95,0.975,0.99),shape=inc_weibull$estimate["shape"],scale=inc_weibull$estimate["scale"]),1)))
par_names <- rbind(names(inc_gamma$estimate), names(inc_lnorm$estimate), names(inc_weibull$estimate))
inc_MS_summary <- cbind(par_names[,1],inc_MS_summary[,1],par_names[,2],inc_MS_summary[,2:ncol(inc_MS_summary)])
write.table(inc_MS_summary,"output/IP_MS_summary.csv",sep=",",col.names = F,row.names = F)

# Save incubation period parameters
inc_params = data.frame(
  inc_ml = coef(inc_lnorm)["meanlog"], inc_sdlog = coef(inc_lnorm)["sdlog"]
)
write.csv(inc_params, "output/inc_params.csv", row.names=FALSE)



## Plot
#-------------

par(mfrow=c(1,1))
par(mar=c(3,3,1,1))

#  Overall incubation period and fit
hist(rabid_carnivores$Incubation.Period.Clean, breaks=seq(-0.5,400.5,5),
     col="lightgrey",main="",xlab="",ylab="",freq=F,axes=F)
axis(1,cex.axis=0.7,padj=-0.5)
axis(2,cex.axis=0.7,padj=0.5)
mtext("Incubation Period",1,2,cex=0.8)
mtext("Density",2,2,cex=0.8)
box(bty="l")
lines(seq(0,400,1),dgamma(seq(0,400,1),shape=coef(inc_gamma)["shape"],rate=coef(inc_gamma)["rate"]),
      col="navy",lwd=1.5)
lines(rep(qgamma(.975,shape=coef(inc_gamma)["shape"],rate=coef(inc_gamma)["rate"]), 100), 
      seq(0,1,length=100), col="navy",lty=1)
lines(seq(0,400,1),dlnorm(seq(0,400,1),meanlog=coef(inc_lnorm)["meanlog"],sdlog=coef(inc_lnorm)["sdlog"]),
      col="orange",lwd=1.5,lty=2)
lines(rep(qlnorm(.975,meanlog=coef(inc_lnorm)["meanlog"],sdlog=coef(inc_lnorm)["sdlog"]), 100), 
      seq(0,1,length=100), col="orange",lty=2)
lines(seq(0,400,1),dweibull(seq(0,400,1),shape=coef(inc_weibull)["shape"],scale=coef(inc_weibull)["scale"]),
      col="green3",lwd=1.5,lty=3)
lines(rep(qweibull(.975,shape=coef(inc_weibull)["shape"],scale=coef(inc_weibull)["scale"]), 100), 
      seq(0,1,length=100), col="green3",lty=3)
lines(rep(quantile(rabid_carnivores$Incubation.Period.Clean,0.975,na.rm=T), 100), 
      seq(0,1,length=100), col="lightgrey",lty=4)
legend("topright",
       paste(c("gamma","lognormal","weibull"),", AIC=",round(c(inc_gamma$aic,inc_lnorm$aic,inc_weibull$aic)),sep=""),
       lty=1:3,col=c("navy","orange","green3"),bty="n",lwd=1.5,cex=0.7)
text(300,0.027,paste("mean=", round(mean(rabid_carnivores$Incubation.Period.Clean,na.rm=T),2),sep=""),cex=0.7)



## Infectious Period
#_______________________

rabid_carnivores$Infectious.Period.Clean <- computePeriod(rabid_carnivores$Infectious.period,
                              rabid_carnivores$Infectious.period.units,
                              rabid_carnivores$Symptoms.started,
                              rep(NA, nrow(rabid_carnivores)), # no end dates possible for infectious periods
                              rabid_carnivores$Symptoms.started.accuracy,
                              rep(NA, nrow(rabid_carnivores)))
summary(rabid_carnivores$Infectious.Period.Clean)

## One case wth infectious period recorded as zero - change to 1 
which(rabid_carnivores$Infectious.Period.Clean==0) # just one case
rabid_carnivores$Infectious.Period.Clean[which(rabid_carnivores$Infectious.Period.Clean==0)]<-1 # change 0 period to 1
table(rabid_carnivores$Infectious.Period.Clean)

## Fit infectious period distributions
inf_idx <- which(!is.na(rabid_carnivores$Infectious.Period.Clean))
inf_gamma = fitdist(rabid_carnivores$Infectious.Period.Clean[inf_idx], "gamma") 
inf_lnorm = fitdist(rabid_carnivores$Infectious.Period.Clean[inf_idx], "lnorm") 
inf_weibull = fitdist(rabid_carnivores$Infectious.Period.Clean[inf_idx], "weibull") 
inf_nbin = fitdist(rabid_carnivores$Infectious.Period.Clean[inf_idx], "nbinom") 
inf_gamma$aic; inf_lnorm$aic; inf_weibull$aic; inf_nbin$aic

## Fit infectious period distributions with censoring to account for use of
## integer days to describe continuous time period
censdata <- data.frame("left"=ifelse(rabid_carnivores$Infectious.Period.Clean[inf_idx]>1,rabid_carnivores$Infectious.Period.Clean[inf_idx]-0.5,0),
                       "right"=ifelse(rabid_carnivores$Infectious.Period.Clean[inf_idx]>1,rabid_carnivores$Infectious.Period.Clean[inf_idx]+0.5,1.5))
inf_cens_gamma <- fitdistcens(censdata, "gamma") 
inf_cens_lnorm <- fitdistcens(censdata, "lnorm") 
inf_cens_weibull <- fitdistcens(censdata, "weibull") 
inf_cens_gamma$aic; inf_cens_lnorm$aic; inf_cens_weibull$aic



## Plot
#-------------

par(mfrow=c(1,2))
par(mar=c(3,3,1,1))

#  Overall Infectious period and fit
hist(rabid_carnivores$Infectious.Period.Clean, breaks=seq(-0.5,15.5,1),cex.main=0.9,
     col="lightgrey",main="Uncensored",xlab="",ylab="",freq=F,axes=F,ylim=c(0,0.36))
axis(1,cex.axis=0.7,padj=-0.5)
axis(2,cex.axis=0.7,padj=0.5)
mtext("Infectious Period",1,2,cex=0.8)
mtext("Density",2,2,cex=0.8)
box(bty="l")
lines(seq(0,20,0.1),dgamma(seq(0,20,0.1),shape=coef(inf_gamma)["shape"],rate=coef(inf_gamma)["rate"]),
      col="navy",lwd=1.5)
lines(rep(qgamma(.975,shape=coef(inf_gamma)["shape"],rate=coef(inf_gamma)["rate"]), 100), 
      seq(0,1,length=100), col="navy",lty=1)
lines(seq(0,20,0.1),dlnorm(seq(0,20,0.1),meanlog=coef(inf_lnorm)["meanlog"],sdlog=coef(inf_lnorm)["sdlog"]),
      col="orange",lwd=1.5,lty=2)
lines(rep(qlnorm(.975,meanlog=coef(inf_lnorm)["meanlog"],sdlog=coef(inf_lnorm)["sdlog"]), 100), 
      seq(0,1,length=100), col="orange",lty=2)
lines(seq(0,20,0.1),dweibull(seq(0,20,0.1),shape=coef(inf_weibull)["shape"],scale=coef(inf_weibull)["scale"]),
      col="green3",lwd=1.5,lty=3)
lines(rep(qweibull(.975,shape=coef(inf_weibull)["shape"],scale=coef(inf_weibull)["scale"]), 100), 
      seq(0,1,length=100), col="green3",lty=3)
points(seq(0,20,1),dnbinom(seq(0,20,1),size=coef(inf_nbin)["size"],mu=coef(inf_nbin)["mu"]),
      col="darkred",lwd=1.5,lty=5,pch=4)
lines(rep(qnbinom(.975,size=coef(inf_nbin)["size"],mu=coef(inf_nbin)["mu"]), 100), 
      seq(0,1,length=100), col="darkred",lty=3)
lines(rep(quantile(rabid_carnivores$Infectious.Period.Clean,0.975,na.rm=T), 100), 
      seq(0,1,length=100), col="lightgrey",lty=4)
legend("topright",
       paste(c("gamma","lognormal","weibull","nbinom"),", AIC=",round(c(inf_gamma$aic,inf_lnorm$aic,inf_weibull$aic)),sep=""),
       lty=c(1:3,NA),pch=c(NA,NA,NA,4),col=c("navy","orange","green3","darkred"),bty="n",lwd=1.5,cex=0.7)
text(10,0.2,paste("mean=", round(mean(rabid_carnivores$Infectious.Period.Clean,na.rm=T),2),sep=""),cex=0.7)

#  Overall Infectious period and fit with interval censoring
hist(rabid_carnivores$Infectious.Period.Clean, breaks=seq(-0.5,15.5,1),cex.main=0.9,
     col="lightgrey",main="Left censored at 1.5",xlab="",ylab="",freq=F,axes=F,ylim=c(0,0.36))
axis(1,cex.axis=0.7,padj=-0.5)
axis(2,cex.axis=0.7,padj=0.5)
mtext("Infectious Period",1,2,cex=0.8)
mtext("Density",2,2,cex=0.8)
box(bty="l")
lines(seq(0,20,0.1),dgamma(seq(0,20,0.1),shape=coef(inf_cens_gamma)["shape"],rate=coef(inf_cens_gamma)["rate"]),
      col="navy",lwd=1.5)
lines(rep(qgamma(.975,shape=coef(inf_cens_gamma)["shape"],rate=coef(inf_cens_gamma)["rate"]), 100), 
      seq(0,1,length=100), col="navy",lty=1)
lines(seq(0,20,0.1),dlnorm(seq(0,20,0.1),meanlog=coef(inf_cens_lnorm)["meanlog"],sdlog=coef(inf_cens_lnorm)["sdlog"]),
      col="orange",lwd=1.5,lty=2)
lines(rep(qlnorm(.975,meanlog=coef(inf_cens_lnorm)["meanlog"],sdlog=coef(inf_cens_lnorm)["sdlog"]), 100), 
      seq(0,1,length=100), col="orange",lty=2)
lines(seq(0,20,0.1),dweibull(seq(0,20,0.1),shape=coef(inf_cens_weibull)["shape"],scale=coef(inf_cens_weibull)["scale"]),
      col="green3",lwd=1.5,lty=3)
lines(rep(qweibull(.975,shape=coef(inf_cens_weibull)["shape"],scale=coef(inf_cens_weibull)["scale"]), 100), 
      seq(0,1,length=100), col="green3",lty=3)
lines(rep(quantile(rabid_carnivores$Infectious.Period.Clean,0.975,na.rm=T), 100), 
      seq(0,1,length=100), col="lightgrey",lty=4)
legend("topright",
       paste(c("gamma","lognormal","weibull"),", AIC=",round(c(inf_cens_gamma$aic,inf_cens_lnorm$aic,inf_cens_weibull$aic)),sep=""),
       lty=1:3,col=c("navy","orange","green3"),bty="n",lwd=1.5,cex=0.7)
text(10,0.2,paste("mean=", round(mean(rabid_carnivores$Infectious.Period.Clean,na.rm=T),2),sep=""),cex=0.7)



## Serial Interval
#_______________________

# Raw serial Intervals
rabid_carnivores$serial_int <- as.numeric(rabid_carnivores$Symptoms.started-rabid_carnivores$Symptoms.started[match(rabid_carnivores$Biter.ID,rabid_carnivores$ID)])
mean(rabid_carnivores$serial_int,na.rm=T)
length(which(!is.na(rabid_carnivores$serial_int)))
summary(rabid_carnivores$serial_int)


# Serial intervals including date uncertainty
rabid_carnivores$Serial.Interval.Clean <- computePeriod(rep(NA, nrow(rabid_carnivores)), # no SIs directly recorded in data
                                                        rep(NA, nrow(rabid_carnivores)), 
                                                        rabid_carnivores$Symptoms.started[match(rabid_carnivores$Biter.ID,rabid_carnivores$ID)],
                                                        rabid_carnivores$Symptoms.started,
                                                        rabid_carnivores$Symptoms.started.uncertainty[match(rabid_carnivores$Biter.ID,rabid_carnivores$ID)],
                                                        rabid_carnivores$Symptoms.started.accuracy)
rabid_carnivores$Serial.Interval.Clean[MultiBite] <- NA # discard intervals for dogs that have multiple possible values
summary(rabid_carnivores$Serial.Interval.Clean)

## Serial intervals simulated from incubation plus infectious period wait time (i.e. time to a bite) 
x = 1000000
incx = rlnorm(x, meanlog = inc_lnorm$estimate["meanlog"], sdlog = inc_lnorm$estimate["sdlog"]) 
infx = rgamma(x, shape = inf_cens_gamma$estimate["shape"], rate = inf_cens_gamma$estimate["rate"])
SI_sims = incx + runif(x, min=0, max=infx)

## 2 cases with serial interval recorded as zero - change to 1 
which(rabid_carnivores$Serial.Interval.Clean==0) 
rabid_carnivores$Serial.Interval.Clean[which(rabid_carnivores$Serial.Interval.Clean==0)]<-1 # change 0 period to 1

## Fit serial interval distributions
ser_idx <- which(!is.na(rabid_carnivores$Serial.Interval.Clean))
ser_gamma = fitdist(rabid_carnivores$Serial.Interval.Clean[ser_idx], "gamma") 
ser_lnorm = fitdist(rabid_carnivores$Serial.Interval.Clean[ser_idx], "lnorm") 
ser_weibull = fitdist(rabid_carnivores$Serial.Interval.Clean[ser_idx], "weibull") 
ser_gamma$aic; ser_lnorm$aic; ser_weibull$aic
boot_ser_gamma <- summary(bootdist(ser_gamma, niter = 1001))
boot_ser_lnorm <- summary(bootdist(ser_lnorm, niter = 1001))
boot_ser_weibull <- summary(bootdist(ser_weibull, niter = 1001))
ser_MS_summary <- rbind(c(paste0(round(ser_gamma$estimate,2)," (",round(boot_ser_gamma$CI[,"2.5%"],2),"-",round(boot_ser_gamma$CI[,"97.5%"],2),")"),round(ser_gamma$aic),round(qgamma(c(0.95,0.975,0.99),shape=ser_gamma$estimate["shape"],rate=ser_gamma$estimate["rate"]),1)),
                        c(paste0(round(ser_lnorm$estimate,2)," (",round(boot_ser_lnorm$CI[,"2.5%"],2),"-",round(boot_ser_lnorm$CI[,"97.5%"],2),")"),round(ser_lnorm$aic),round(qlnorm(c(0.95,0.975,0.99),meanlog=ser_lnorm$estimate["meanlog"],sdlog=ser_lnorm$estimate["sdlog"]),1)),
                        c(paste0(round(ser_weibull$estimate,2)," (",round(boot_ser_weibull$CI[,"2.5%"],2),"-",round(boot_ser_weibull$CI[,"97.5%"],2),")"),round(ser_weibull$aic),round(qweibull(c(0.95,0.975,0.99),shape=ser_weibull$estimate["shape"],scale=ser_weibull$estimate["scale"]),1)))
par_names <- rbind(names(ser_gamma$estimate), names(ser_lnorm$estimate), names(ser_weibull$estimate))
ser_MS_summary <- cbind(par_names[,1],ser_MS_summary[,1],par_names[,2],ser_MS_summary[,2:ncol(ser_MS_summary)])
write.table(ser_MS_summary,"output/SI_MS_summary.csv",sep=",",col.names = F,row.names = F)

## Simulate serial interval convolutions
SIdouble_lnorm_sim <- rlnorm(x, meanlog = ser_lnorm$estimate["meanlog"], sdlog = ser_lnorm$estimate["sdlog"]) + rlnorm(x, meanlog = ser_lnorm$estimate["meanlog"], sdlog = ser_lnorm$estimate["sdlog"])
SIdouble_lnorm <- fitdist(SIdouble_lnorm_sim,"lnorm")
SIdouble_gamma <- fitdist(SIdouble_lnorm_sim,"gamma")
SIdouble_weibull <- fitdist(SIdouble_lnorm_sim,"weibull")
SIdouble_lnorm$aic; SIdouble_gamma$aic;SIdouble_weibull$aic # lognormal works best here

# Write parameters for transmission tree inference
SI_params = data.frame(
  SI_ml = coef(ser_lnorm)["meanlog"], SI_sdlog = coef(ser_lnorm)["sdlog"],
  SI2_ml = coef(SIdouble_lnorm)["meanlog"], SI2_sdlog = coef(SIdouble_lnorm)["sdlog"]
)
write.csv(SI_params, "output/SI_params.csv", row.names=FALSE)

# Write serial intervals for plotting
write.table(rabid_carnivores$Serial.Interval.Clean[!is.na(rabid_carnivores$Serial.Interval.Clean)],
          "output/serial_intervals.csv", row.names = F, col.names=F)


## Plots
#-------------

# pdf("figs/SerialInterval.pdf",7,6)
par(mfrow=c(2,2))
par(mar=c(3,3,2,1))

# Raw serial intervals
hist(rabid_carnivores$serial_int, breaks=seq(-0.5,420.5,5),cex.main=0.9,
     col="lightgrey",main="Raw serial intervals",xlab="",ylab="",freq=F,axes=F)
axis(1,cex.axis=0.7,padj=-0.5)
axis(2,cex.axis=0.7,padj=0.5)
mtext("Serial Interval",1,2,cex=0.8)
mtext("Density",2,2,cex=0.8)
box(bty="l")

#  Cleaned serial intervals and fit
hist(rabid_carnivores$Serial.Interval.Clean, breaks=seq(-0.5,420.5,5),cex.main=0.9,
     col="lightgrey",main="Cleaned serial intervals",xlab="",ylab="",freq=F,axes=F,
     ylim=c(0,0.04))
axis(1,cex.axis=0.7,padj=-0.5)
axis(2,cex.axis=0.7,padj=0.5)
mtext("Serial Interval",1,2,cex=0.8)
mtext("Density",2,2,cex=0.8)
box(bty="l")
lines(seq(0,450,1),dgamma(seq(0,450,1),shape=coef(ser_gamma)["shape"],rate=coef(ser_gamma)["rate"]),
      col="navy",lwd=1.5)
lines(rep(qgamma(.975,shape=coef(ser_gamma)["shape"],rate=coef(ser_gamma)["rate"]), 100), 
      seq(0,1,length=100), col="navy",lty=1)
lines(seq(0,450,1),dlnorm(seq(0,450,1),meanlog=coef(ser_lnorm)["meanlog"],sdlog=coef(ser_lnorm)["sdlog"]),
      col="orange",lwd=1.5,lty=2)
lines(rep(qlnorm(.975,meanlog=coef(ser_lnorm)["meanlog"],sdlog=coef(ser_lnorm)["sdlog"]), 100), 
      seq(0,1,length=100), col="orange",lty=2)
lines(seq(0,450,1),dweibull(seq(0,450,1),shape=coef(ser_weibull)["shape"],scale=coef(ser_weibull)["scale"]),
      col="green3",lwd=1.5,lty=3)
lines(rep(qweibull(.975,shape=coef(ser_weibull)["shape"],scale=coef(ser_weibull)["scale"]), 100), 
      seq(0,1,length=100), col="green3",lty=3)
lines(rep(quantile(rabid_carnivores$Serial.Interval.Clean,0.975,na.rm=T), 100), 
      seq(0,1,length=100), col="lightgrey",lty=4)
legend("topright",
       paste(c("gamma","lognormal","weibull"),", AIC=",round(c(ser_gamma$aic,ser_lnorm$aic,ser_weibull$aic)),sep=""),
       lty=c(1:3),col=c("navy","orange","green3"),bty="n",lwd=1.5,cex=0.7)
text(350,0.03,paste("mean=", round(mean(rabid_carnivores$Serial.Interval.Clean,na.rm=T),2),sep=""),cex=0.7)

# Serial intervals simulated from incubation period & time to bite
hist(SI_sims, breaks=seq(-0.5,3015.5,5),cex.main=0.9,xlim=c(0,420),ylim=c(0,0.04),
     col="lightgrey",main="Serial intervals simulated from incubation period\n& time to bite",xlab="",ylab="",freq=F,axes=F)
axis(1,cex.axis=0.7,padj=-0.5)
axis(2,cex.axis=0.7,padj=0.5)
mtext("Serial Interval",1,2,cex=0.8)
mtext("Density",2,2,cex=0.8)
box(bty="l")
lines(seq(0,450,1),dlnorm(seq(0,450,1),meanlog=coef(ser_lnorm)["meanlog"],sdlog=coef(ser_lnorm)["sdlog"]),
      col="orange",lwd=1.5,lty=2)
text(350,0.03,paste("mean=", round(mean(SI_sims),2),sep=""),cex=0.7)



#  Double serial interval
hist(SIdouble_lnorm_sim, breaks=seq(-0.5,3300.5,5),cex.main=0.9,xlim=c(0,420),ylim=c(0,0.02),
     col="lightgrey",main="Double SI (with lnorm fit)",xlab="",ylab="",freq=F,axes=F)
axis(1,cex.axis=0.7,padj=-0.5)
axis(2,cex.axis=0.7,padj=0.5)
mtext("Double Serial Interval",1,2,cex=0.8)
mtext("Density",2,2,cex=0.8)
box(bty="l")
lines(seq(0,450,1),dlnorm(seq(0,450,1),meanlog=coef(SIdouble_lnorm)["meanlog"],sdlog=coef(SIdouble_lnorm)["sdlog"]),
      col="orange",lwd=1.5,lty=2)
lines(rep(qlnorm(.975,meanlog=coef(SIdouble_lnorm)["meanlog"],sdlog=coef(SIdouble_lnorm)["sdlog"]), 100), 
      seq(0,1,length=100), col="orange",lty=2)
lines(rep(quantile(SIdouble_lnorm_sim,0.975,na.rm=T), 100), 
      seq(0,1,length=100), col="lightgrey",lty=4)

text(350,0.015,paste("mean=", round(mean(SIdouble_lnorm_sim),2),sep=""),cex=0.7)

# dev.off()



## Generation Interval
#_______________________

# Raw generation Interval
rabid_carnivores$generation_int <- as.numeric(rabid_carnivores$Date.bitten-rabid_carnivores$Date.bitten[match(rabid_carnivores$Biter.ID,rabid_carnivores$ID)])
mean(rabid_carnivores$generation_int,na.rm=T)
length(which(!is.na(rabid_carnivores$generation_int)))


# Generation intervals including date uncertainty
rabid_carnivores$Generation.Interval.Clean <- computePeriod(rep(NA, nrow(rabid_carnivores)), # no GIs directly recorded in data
                                                            rep(NA, nrow(rabid_carnivores)), 
                                                            rabid_carnivores$Date.bitten[match(rabid_carnivores$Biter.ID,rabid_carnivores$ID)],
                                                            rabid_carnivores$Date.bitten,
                                                            rabid_carnivores$Date.bitten.uncertainty[match(rabid_carnivores$Biter.ID,rabid_carnivores$ID)],
                                                            rabid_carnivores$Date.bitten.uncertainty)
summary(rabid_carnivores$Generation.Interval.Clean)


## Fit generation interval distributions
gen_idx <- which(!is.na(rabid_carnivores$Generation.Interval.Clean))
gen_gamma = fitdist(rabid_carnivores$Generation.Interval.Clean[gen_idx], "gamma") 
gen_lnorm = fitdist(rabid_carnivores$Generation.Interval.Clean[gen_idx], "lnorm") 
gen_weibull = fitdist(rabid_carnivores$Generation.Interval.Clean[gen_idx], "weibull") 
gen_gamma$aic; gen_lnorm$aic; gen_weibull$aic


# Generation intervals taking only shortest one for each biter
rabid_carnivores$smallestGI <- NA
length(which(rabid_carnivores$Biter.ID==0 & !is.na(rabid_carnivores$Date.bitten))) #43 of these where biter unknown but none will have GIs
for(i in unique(rabid_carnivores$Biter.ID)){
  biter_GIs_available <- which(rabid_carnivores$Biter.ID==i & !is.na(rabid_carnivores$generation_int))
  if(length(biter_GIs_available>0) & i!=0){
    rabid_carnivores$smallestGI[biter_GIs_available[which.min(rabid_carnivores$generation_int[biter_GIs_available])]] <- 1
  }else if(i==0){
    rabid_carnivores$smallestGI[biter_GIs_available] <- 1
  }
}
summary(rabid_carnivores$Generation.Interval.Clean[which(rabid_carnivores$smallestGI==1)])


## Fit shortest generation interval distributions
gen_short_idx <- which(!is.na(rabid_carnivores$Generation.Interval.Clean) & rabid_carnivores$smallestGI==1)
gen_short_gamma = fitdist(rabid_carnivores$Generation.Interval.Clean[gen_short_idx], "gamma") 
gen_short_lnorm = fitdist(rabid_carnivores$Generation.Interval.Clean[gen_short_idx], "lnorm") 
gen_short_weibull = fitdist(rabid_carnivores$Generation.Interval.Clean[gen_short_idx], "weibull") 
gen_short_gamma$aic; gen_short_lnorm$aic; gen_short_weibull$aic



## Plot
#-------------

# pdf("figs/GenerationInterval.pdf",7,6)
par(mfrow=c(2,2))
par(mar=c(3,3,2,1))

# Raw generation intervals
hist(rabid_carnivores$generation_int, breaks=seq(-0.5,250.5,5),cex.main=0.9,
     col="lightgrey",main="Raw generation intervals",xlab="",ylab="",freq=F,axes=F)
axis(1,cex.axis=0.7,padj=-0.5)
axis(2,cex.axis=0.7,padj=0.5)
mtext("Generation Interval",1,2,cex=0.8)
mtext("Density",2,2,cex=0.8)
box(bty="l")
text(200,0.025,paste("mean=", round(mean(rabid_carnivores$generation_int,na.rm=T),2),sep=""),cex=0.7)


# Cleaned generation intervals and fit
hist(rabid_carnivores$Generation.Interval.Clean, breaks=seq(-0.5,250.5,5),cex.main=0.9,
     col="lightgrey",main="Cleaned generation intervals",xlab="",ylab="",freq=F,axes=F,
     ylim=c(0,0.035))
axis(1,cex.axis=0.7,padj=-0.5)
axis(2,cex.axis=0.7,padj=0.5)
mtext("Generation Interval",1,2,cex=0.8)
mtext("Density",2,2,cex=0.8)
box(bty="l")
lines(seq(0,450,1),dgamma(seq(0,450,1),shape=coef(gen_gamma)["shape"],rate=coef(gen_gamma)["rate"]),
      col="navy",lwd=1.5)
lines(rep(qgamma(.975,shape=coef(gen_gamma)["shape"],rate=coef(gen_gamma)["rate"]), 100), 
      seq(0,1,length=100), col="navy",lty=1)
lines(seq(0,450,1),dlnorm(seq(0,450,1),meanlog=coef(gen_lnorm)["meanlog"],sdlog=coef(gen_lnorm)["sdlog"]),
      col="orange",lwd=1.5,lty=2)
lines(rep(qlnorm(.975,meanlog=coef(gen_lnorm)["meanlog"],sdlog=coef(gen_lnorm)["sdlog"]), 100), 
      seq(0,1,length=100), col="orange",lty=2)
lines(seq(0,450,1),dweibull(seq(0,450,1),shape=coef(gen_weibull)["shape"],scale=coef(gen_weibull)["scale"]),
      col="green3",lwd=1.5,lty=3)
lines(rep(qweibull(.975,shape=coef(gen_weibull)["shape"],scale=coef(gen_weibull)["scale"]), 100), 
      seq(0,1,length=100), col="green3",lty=3)
lines(rep(quantile(rabid_carnivores$Generation.Interval.Clean,0.975,na.rm=T), 100), 
      seq(0,1,length=100), col="lightgrey",lty=4)
legend("topright",
       paste(c("gamma","lognormal","weibull"),", AIC=",round(c(gen_gamma$aic,gen_lnorm$aic,gen_weibull$aic)),sep=""),
       lty=c(1:3),col=c("navy","orange","green3"),bty="n",lwd=1.5,cex=0.7)
text(200,0.025,paste("mean=", round(mean(rabid_carnivores$Generation.Interval.Clean,na.rm=T),2),sep=""),cex=0.7)


# Cleaned first generation intervals and fit
hist(rabid_carnivores$Generation.Interval.Clean[which(rabid_carnivores$smallestGI==1)], breaks=seq(-0.5,250.5,5),cex.main=0.9,
     col="lightgrey",main="Cleaned first generation interval for each dog",xlab="",ylab="",freq=F,axes=F,
     ylim=c(0,0.035))
axis(1,cex.axis=0.7,padj=-0.5)
axis(2,cex.axis=0.7,padj=0.5)
mtext("Generation Interval",1,2,cex=0.8)
mtext("Density",2,2,cex=0.8)
box(bty="l")
lines(seq(0,450,1),dgamma(seq(0,450,1),shape=coef(gen_short_gamma)["shape"],rate=coef(gen_short_gamma)["rate"]),
      col="navy",lwd=1.5)
lines(rep(qgamma(.975,shape=coef(gen_short_gamma)["shape"],rate=coef(gen_short_gamma)["rate"]), 100), 
      seq(0,1,length=100), col="navy",lty=1)
lines(seq(0,450,1),dlnorm(seq(0,450,1),meanlog=coef(gen_short_lnorm)["meanlog"],sdlog=coef(gen_short_lnorm)["sdlog"]),
      col="orange",lwd=1.5,lty=2)
lines(rep(qlnorm(.975,meanlog=coef(gen_short_lnorm)["meanlog"],sdlog=coef(gen_short_lnorm)["sdlog"]), 100), 
      seq(0,1,length=100), col="orange",lty=2)
lines(seq(0,450,1),dweibull(seq(0,450,1),shape=coef(gen_short_weibull)["shape"],scale=coef(gen_short_weibull)["scale"]),
      col="green3",lwd=1.5,lty=3)
lines(rep(qweibull(.975,shape=coef(gen_short_weibull)["shape"],scale=coef(gen_short_weibull)["scale"]), 100), 
      seq(0,1,length=100), col="green3",lty=3)
lines(rep(quantile(rabid_carnivores$Generation.Interval.Clean,0.975,na.rm=T), 100), 
      seq(0,1,length=100), col="lightgrey",lty=4)
legend("topright",
       paste(c("gamma","lognormal","weibull"),", AIC=",round(c(gen_short_gamma$aic,gen_short_lnorm$aic,gen_short_weibull$aic)),sep=""),
       lty=c(1:3),col=c("navy","orange","green3"),bty="n",lwd=1.5,cex=0.7)
text(200,0.025,paste("mean=", round(mean(rabid_carnivores$Generation.Interval.Clean[which(rabid_carnivores$smallestGI==1)],na.rm=T),2),sep=""),cex=0.7)

# dev.off()





