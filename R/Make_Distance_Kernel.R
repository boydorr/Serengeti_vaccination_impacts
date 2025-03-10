rm(list=ls())
library(fitdistrplus)

set.seed(0)



# Read in data on rabid carnivores
rabid_carnivores <- readRDS(file = "output/clean_bite_data.rda")
nrow(rabid_carnivores) 

# Read in all Serengeti CT data
ct <- readRDS("data/animalCT_deid.rds")

# Subset on location information being available
rabid_carnivores_location <- subset(rabid_carnivores, !is.na(UTM.Easting))
nrow(rabid_carnivores_location) 

## calculate max distance moved by biter (locations of bitten & origin)
max_moved=function(biter, bitten){
  gps=rbind(cbind(biter$UTM.Easting, biter$UTM.Northing),
            cbind(bitten$UTM.Easting, bitten$UTM.Northing))
  max(as.matrix(dist(gps)),na.rm=T)
}

## all distances moved by biter
all_moved=function(biter, bitten){
  gps=rbind(cbind(biter$UTM.Easting, biter$UTM.Northing),
            cbind(bitten$UTM.Easting, bitten$UTM.Northing))
  distances <- as.numeric(as.matrix(dist(gps))[1,2:(nrow(bitten)+1)])
  distances[which(!is.na(distances))]
}



## calculate case-by-case and determine whether needs to be censored or not
rabid_carnivores$max_dist_case = rabid_carnivores$max_dist_contact = NA
rabid_carnivores$max_dist_case_censored = rabid_carnivores$max_dist_contact_censored = NA
dist_case_all <- dist_contact_all <- data.frame("dist"=numeric(),"right_censored"=logical(),"ID"=integer())
for(i in 1:nrow(rabid_carnivores)){
  cases = which(rabid_carnivores$Biter.ID == rabid_carnivores$ID[i])
  contacts = which(ct$Biter.ID == rabid_carnivores$ID[i])
  rabid_carnivores$max_dist_case[i] = ifelse(length(cases)>0, max_moved(rabid_carnivores[i,], rabid_carnivores[cases,]), NA)
  rabid_carnivores$max_dist_contact[i] = ifelse(length(contacts)>0, max_moved(rabid_carnivores[i,], ct[contacts,]), NA)
  if(rabid_carnivores$Owner[i]!="Known" & !is.na(rabid_carnivores$max_dist_case[i])){rabid_carnivores$max_dist_case_censored[i]<-T
  }else if(rabid_carnivores$Owner[i]=="Known" & !is.na(rabid_carnivores$max_dist_case[i])){rabid_carnivores$max_dist_case_censored[i]<-F}
  if(rabid_carnivores$Owner[i]!="Known" & !is.na(rabid_carnivores$max_dist_contact[i])){rabid_carnivores$max_dist_contact_censored[i]<-T
  }else if(rabid_carnivores$Owner[i]=="Known" & !is.na(rabid_carnivores$max_dist_contact[i])){rabid_carnivores$max_dist_contact_censored[i]<-F}
  if(length(cases)>0){
    all_case_moves <- all_moved(rabid_carnivores[i,], rabid_carnivores[cases,])
    if(rabid_carnivores$Owner[i]!="Known"){right_censored <- rep(T,length(all_case_moves))
    }else{right_censored <- rep(F,length(all_case_moves))}
    dist_case_all <- rbind(dist_case_all,cbind(all_case_moves,right_censored,rep(rabid_carnivores$ID[i],length(all_case_moves))))
  }
  if(length(contacts)>0){
    all_contact_moves <- all_moved(rabid_carnivores[i,], ct[contacts,])
    if(rabid_carnivores$Owner[i]!="Known"){right_censored <- rep(T,length(all_contact_moves))
    }else{right_censored <- rep(F,length(all_contact_moves))}
    dist_contact_all <- rbind(dist_contact_all,cbind(all_contact_moves,right_censored,rep(rabid_carnivores$ID[i],length(all_contact_moves))))
  }
  # print(i)
}
names(dist_contact_all)<-names(dist_case_all)<-c("distance","right_censored","ID")
dist_contact_all$distance <- as.numeric(dist_contact_all$distance)
dist_contact_all$right_censored <- as.logical(dist_contact_all$right_censored)
dist_contact_all$ID <- as.integer(dist_contact_all$distance)
dist_case_all$distance <- as.numeric(dist_case_all$distance)
dist_case_all$right_censored <- as.logical(dist_case_all$right_censored)
dist_case_all$ID <- as.integer(dist_case_all$distance)


# Compare distances between cases and contacts
par(mfrow=c(1,2)) 
hist(rabid_carnivores$max_dist_case, breaks=seq(0, 25000, 100), xlim=c(0,5000),freq=F)
hist(rabid_carnivores$max_dist_contact, breaks=seq(0, 25000, 100), xlim=c(0,5000),freq=F) 
mean(rabid_carnivores$max_dist_case, na.rm = TRUE); median(rabid_carnivores$max_dist_case, na.rm = TRUE) 
mean(rabid_carnivores$max_dist_contact, na.rm = TRUE); median(rabid_carnivores$max_dist_contact, na.rm = TRUE) 
hist(dist_case_all$distance, breaks=seq(0, 25000, 100), xlim=c(0,5000),freq=F)
hist(dist_contact_all$distance, breaks=seq(0, 25000, 100), xlim=c(0,5000),freq=F) 
mean(dist_case_all$distance, na.rm = TRUE); median(dist_case_all$distance, na.rm = TRUE) 
mean(dist_contact_all$distance, na.rm = TRUE);  median(dist_contact_all$distance, na.rm = TRUE) 

# Compare distances between contacts for unknown and known biters
mean(rabid_carnivores$max_dist_contact, na.rm = TRUE); quantile(rabid_carnivores$max_dist_contact,c(0.5,0.95,0.975), na.rm = TRUE) 
mean(rabid_carnivores$max_dist_contact[which(rabid_carnivores$max_dist_contact_censored==F)], na.rm = TRUE);  quantile(rabid_carnivores$max_dist_contact[which(rabid_carnivores$max_dist_contact_censored==F)],c(0.5,0.95,0.975), na.rm = TRUE) 
mean(rabid_carnivores$max_dist_contact[which(rabid_carnivores$max_dist_contact_censored==F|rabid_carnivores$max_dist_contact>100)], na.rm = TRUE);  quantile(rabid_carnivores$max_dist_contact[which(rabid_carnivores$max_dist_contact_censored==F|rabid_carnivores$max_dist_contact>100)],c(0.5,0.95,0.975), na.rm = TRUE) 
mean(ifelse(rabid_carnivores$Owner=="Known",rabid_carnivores$max_dist_contact,rabid_carnivores$max_dist_contact+100), na.rm = TRUE);  quantile(ifelse(rabid_carnivores$Owner=="Known",rabid_carnivores$max_dist_contact,rabid_carnivores$max_dist_contact+100),c(0.5,0.95,0.975), na.rm = TRUE) 

mean(dist_contact_all$distance, na.rm = TRUE); quantile(dist_contact_all$distance,c(0.5,0.95,0.975), na.rm = TRUE) 
mean(dist_contact_all$distance[which(dist_contact_all$right_censored==F)], na.rm = TRUE);  quantile(dist_contact_all$distance[which(dist_contact_all$right_censored==F)],c(0.5,0.95,0.975), na.rm = TRUE) 
mean(dist_contact_all$distance[which(dist_contact_all$right_censored==F|dist_contact_all$distance>100)], na.rm = TRUE);  quantile(dist_contact_all$distance[which(dist_contact_all$right_censored==F|dist_contact_all$distance>100)],c(0.5,0.95,0.975), na.rm = TRUE) 
mean(ifelse(dist_contact_all$right_censored==T,dist_contact_all$distance,dist_contact_all$distance+100), na.rm = TRUE);  quantile(ifelse(dist_contact_all$right_censored==T,dist_contact_all$distance,dist_contact_all$distance+100),c(0.5,0.95,0.975), na.rm = TRUE) 
mean(dist_contact_all$distance[which(dist_contact_all$right_censored==F|(dist_contact_all$ID %in% unique(dist_contact_all$ID[which(dist_contact_all$distance>100)])))],na.rm=T); quantile(dist_contact_all$distance[which(dist_contact_all$right_censored==F|(dist_contact_all$ID %in% unique(dist_contact_all$ID[which(dist_contact_all$distance>100)])))],c(0.5,0.95,0.975),na.rm=T)


# Fit distributions to maximum distances
idx<-which(!is.na(rabid_carnivores$max_dist_contact))
censdata <- data.frame("left"=ifelse(rabid_carnivores$max_dist_contact[idx]>100,rabid_carnivores$max_dist_contact[idx],0),
                       "right"=ifelse(rabid_carnivores$max_dist_contact[idx]>100,rabid_carnivores$max_dist_contact[idx],100))
censdata$left[which(rabid_carnivores$max_dist_contact_censored[idx]==T)] <- ifelse(censdata$left[which(rabid_carnivores$max_dist_contact_censored[idx]==T)]<100,100,censdata$left[which(rabid_carnivores$max_dist_contact_censored[idx]==T)])
censdata$right[which(rabid_carnivores$max_dist_contact_censored[idx]==T)] <- NA
max_dist_gamma <- fitdistcens(censdata, "gamma",control = list(ndeps=c(1e-4,1e-4), parscale=c(1,1e-5))) 
max_dist_lnorm <- fitdistcens(censdata, "lnorm") 
max_dist_weibull <- fitdistcens(censdata, "weibull") 
max_dist_gamma$aic; max_dist_lnorm$aic; max_dist_weibull$aic # weibull best (but only by AIC diff ~1)
qweibull(c(0.95,0.975,0.99),shape=max_dist_weibull$estimate["shape"],scale=max_dist_weibull$estimate["scale"])


# Fit distributions to all distances
censdata <- data.frame("left"=ifelse(dist_contact_all$distance>100,dist_contact_all$distance,0),
                       "right"=ifelse(dist_contact_all$distance>100,dist_contact_all$distance,100))
censdata$left[which(dist_contact_all$right_censored==T)] <- ifelse(censdata$left[which(dist_contact_all$right_censored==T)]<100,100,censdata$left[which(dist_contact_all$right_censored==T)])
censdata$right[which(dist_contact_all$right_censored==T)] <- NA
all_dist_gamma <- fitdistcens(censdata, "gamma",control = list(ndeps=c(1e-4,1e-4), parscale=c(1,1e-5))) 
all_dist_lnorm <- fitdistcens(censdata, "lnorm") 
all_dist_weibull <- fitdistcens(censdata, "weibull") 
all_dist_gamma$aic; all_dist_lnorm$aic; all_dist_weibull$aic # weibull best
boot_all_dist_gamma <- summary(bootdistcens(all_dist_gamma, niter = 1001))
boot_all_dist_lnorm <- summary(bootdistcens(all_dist_lnorm, niter = 1001))
boot_all_dist_weibull <- summary(bootdistcens(all_dist_weibull, niter = 1001))
all_dist_MS_summary <- rbind(c(paste0(signif(all_dist_gamma$estimate,2)," (",signif(boot_all_dist_gamma$CI[,"2.5%"],2),"-",signif(boot_all_dist_gamma$CI[,"97.5%"],2),")"),round(all_dist_gamma$aic),round(qgamma(c(0.95,0.975,0.99),shape=all_dist_gamma$estimate["shape"],rate=all_dist_gamma$estimate["rate"]),1)),
                        c(paste0(round(all_dist_lnorm$estimate,2)," (",round(boot_all_dist_lnorm$CI[,"2.5%"],2),"-",round(boot_all_dist_lnorm$CI[,"97.5%"],2),")"),round(all_dist_lnorm$aic),round(qlnorm(c(0.95,0.975,0.99),meanlog=all_dist_lnorm$estimate["meanlog"],sdlog=all_dist_lnorm$estimate["sdlog"]),1)),
                        c(paste0(round(all_dist_weibull$estimate,2)," (",round(boot_all_dist_weibull$CI[,"2.5%"],2),"-",round(boot_all_dist_weibull$CI[,"97.5%"],2),")"),round(all_dist_weibull$aic),round(qweibull(c(0.95,0.975,0.99),shape=all_dist_weibull$estimate["shape"],scale=all_dist_weibull$estimate["scale"]),1)))
par_names <- rbind(names(all_dist_gamma$estimate), names(all_dist_lnorm$estimate), names(all_dist_weibull$estimate))
all_dist_MS_summary <- cbind(par_names[,1],all_dist_MS_summary[,1],par_names[,2],all_dist_MS_summary[,2:ncol(all_dist_MS_summary)])
# write.table(all_dist_MS_summary,"output/DK_MS_summary.csv",sep=",",col.names = F,row.names = F)



# Convolution based on max distances (using base distribution with lowest AIC)
n = 1000000
sim_max_dists_2 <- rweibull(n,shape=coef(max_dist_weibull)["shape"],scale=coef(max_dist_weibull)["scale"])+rweibull(n,shape=coef(max_dist_weibull)["shape"],scale=coef(max_dist_weibull)["scale"])
max_dist_gamma_2 <- fitdist(sim_max_dists_2,"gamma",control = list(ndeps=c(1e-4,1e-4), parscale=c(1,1e-5)))
max_dist_weibull_2 <- fitdist(sim_max_dists_2,"weibull")
max_dist_lnorm_2 <- fitdist(sim_max_dists_2,"lnorm")
max_dist_gamma_2$aic; max_dist_lnorm_2$aic; max_dist_weibull_2$aic # weibull best

# Convolution based on all distances (using base distribution with lowest AIC)
sim_all_dists_2 <- rweibull(n,shape=coef(all_dist_weibull)["shape"],scale=coef(all_dist_weibull)["scale"])+rweibull(n,shape=coef(all_dist_weibull)["shape"],scale=coef(all_dist_weibull)["scale"])
all_dist_gamma_2 <- fitdist(sim_all_dists_2,"gamma",control = list(ndeps=c(1e-4,1e-4), parscale=c(1,1e-5)))
all_dist_weibull_2 <- fitdist(sim_all_dists_2,"weibull")
all_dist_lnorm_2 <- fitdist(sim_all_dists_2,"lnorm")
all_dist_gamma_2$aic; all_dist_lnorm_2$aic; all_dist_weibull_2$aic # gamma best


# But this assumes that the subsequent distances are in the same direction -
# probably makes more sense to assume a random change in direction between
# distances

# Simulate two max distances with an intermediate random change in direction
conv_with_direction <- function(n=1000000,dist,par1,par2){
  
  angle1 <- runif(n = n, min = 0, max = 2*pi)  
  angle2 <- runif(n = n, min = 0, max = 2*pi)  
  dist_1 <- get(paste0("r",dist))(n,par1,par2)
  dist_2 <- get(paste0("r",dist))(n,par1,par2)
  sim_2_dists <- rep(NA,n)
  
  for(i in 1:n){
    
    # Coords after first distance
    x1 <- (sin(angle1[i]) * dist_1[i]) + 0  
    y1 <- (cos(angle1[i]) * dist_1[i]) + 0
    
    # Coords after second distance
    x2 <- (sin(angle2[i]) * dist_2[i]) + x1  
    y2 <- (cos(angle2[i]) * dist_2[i]) + y1
    
    # Total distance travelled from (0,0)
    sim_2_dists[i] <- sqrt(x2^2 + y2^2)
    
  }
  
  return(sim_2_dists)
}
sim_max_dists_2_direction <- conv_with_direction(dist="weibull",par1=coef(max_dist_weibull)["shape"],par2=coef(max_dist_weibull)["scale"])
summary(sim_max_dists_2_direction)
max_dist_gamma_2_direction <- fitdist(sim_max_dists_2_direction,"gamma",control = list(ndeps=c(1e-4,1e-4), parscale=c(1,1e-5)))
max_dist_lnorm_2_direction <- fitdist(sim_max_dists_2_direction,"lnorm")
max_dist_weibull_2_direction <- fitdist(sim_max_dists_2_direction,"weibull")
max_dist_gamma_2_direction$aic; max_dist_lnorm_2_direction$aic; max_dist_weibull_2_direction$aic # weibull best


# Simulate two draws from weibull with an intermediate random change in direction in between
sim_dists_2_direction <- conv_with_direction(dist="weibull",par1=coef(all_dist_weibull)["shape"],par2=coef(all_dist_weibull)["scale"])
summary(sim_dists_2_direction)
all_dist_gamma_2_direction <- fitdist(sim_dists_2_direction,"gamma",control = list(ndeps=c(1e-4,1e-4), parscale=c(1,1e-5)))
all_dist_lnorm_2_direction <- fitdist(sim_dists_2_direction,"lnorm")
all_dist_weibull_2_direction <- fitdist(sim_dists_2_direction,"weibull")
all_dist_gamma_2_direction$aic; all_dist_lnorm_2_direction$aic; all_dist_weibull_2_direction$aic # weibull best

# Simulate two draws from gamma with an intermediate random change in direction in between
sim_dists_2_direction_gamma <- conv_with_direction(dist="gamma",par1=coef(all_dist_gamma)["shape"],par2=coef(all_dist_gamma)["rate"])
summary(sim_dists_2_direction_gamma)
all_dist_gamma_2_direction_gamma <- fitdist(sim_dists_2_direction_gamma,"gamma",control = list(ndeps=c(1e-4,1e-4), parscale=c(1,1e-5)))
all_dist_lnorm_2_direction_gamma <- fitdist(sim_dists_2_direction_gamma,"lnorm")
all_dist_weibull_2_direction_gamma <- fitdist(sim_dists_2_direction_gamma,"weibull")
all_dist_gamma_2_direction_gamma$aic; all_dist_lnorm_2_direction_gamma$aic; all_dist_weibull_2_direction_gamma$aic # gamma best

# Simulate two draws from lnorm with an intermediate random change in direction in between
sim_dists_2_direction_lnorm <- conv_with_direction(dist="lnorm",par1=coef(all_dist_lnorm)["meanlog"],par2=coef(all_dist_lnorm)["sdlog"])
summary(sim_dists_2_direction_lnorm)
all_dist_gamma_2_direction_lnorm <- fitdist(sim_dists_2_direction_lnorm,"gamma",control = list(ndeps=c(1e-4,1e-4), parscale=c(1,1e-5)))
all_dist_lnorm_2_direction_lnorm <- fitdist(sim_dists_2_direction_lnorm,"lnorm")
all_dist_weibull_2_direction_lnorm <- fitdist(sim_dists_2_direction_lnorm,"weibull")
all_dist_gamma_2_direction_lnorm$aic; all_dist_lnorm_2_direction_lnorm$aic; all_dist_weibull_2_direction_lnorm$aic # lnorm best


# Write parameters for transmission tree inference
# (This section is commented out to avoid overwriting the distributions and 
# distances actually used in the paper, which used the exact case coordinates,
# not the ones above that were deidentified by jittering)
# summarise_dist <- function(dist_fit){
#   
#   c(dist = dist_fit$distname,
#     AIC = round(dist_fit$aic),
#     v1name = names(dist_fit$estimate)[1], 
#     var1 = dist_fit$estimate[1],
#     v2name = names(dist_fit$estimate)[2],
#     var2 = dist_fit$estimate[2])
#   
# }
# dists <- list(all_dist_gamma, all_dist_gamma_2_direction_gamma, 
#               all_dist_lnorm, all_dist_lnorm_2_direction_lnorm,
#               all_dist_weibull, all_dist_weibull_2_direction)
# DK_params = data.frame(
#   kernel_type=rep(c("DK","DK2"),3),
#   dist=NA,
#   AIC=NA,
#   par1name=NA,
#   par1est=NA,
#   par2name=NA,
#   par2est=NA
# )
# for(i in 1:length(dists)){
#     DK_params[i,2:ncol(DK_params)] <- summarise_dist(dists[[i]])
# }
# DK_params$AIC[which(DK_params$kernel_type=="DK2")]<-NA
# write.csv(DK_params, "output/DK_params.csv", row.names=FALSE)
# 
# 
# # Write distances fitted to for plotting
# write.table(dist_contact_all$distance,
#             "output/distances_between_case_contacts.csv", row.names = F, col.names=F)




## Plot
#-----------

hist(dist_contact_all$distance, breaks=seq(-0.5,35000.5,200),cex.main=0.9,
     col="dodgerblue",border="dodgerblue",main="",xlab="",ylab="",freq=F,axes=F,xlim=c(0,10000))
axis(1,cex.axis=0.9,padj=-0.5)
axis(2,cex.axis=0.9,padj=0.5)
mtext("Distance from biter to bitee",1,2,cex=1.4)
mtext("Density",2,2,cex=1.4)
box(bty="l")
lines(rep(quantile(dist_contact_all$distance,0.975,na.rm=T), 100), 
      seq(0,1,length=100), col="red",lty=4,lwd=2)



# pdf("figs/SpatialInfectionKernel.pdf",7,9)
# par(mfrow=c(3,2))
# par(mar=c(3,3,3,1))
# 
# #  Kernel fitted to max distances
# hist(rabid_carnivores$max_dist_contact, breaks=seq(-0.5,35000.5,500),cex.main=0.9,
#      col="lightgrey",main="Using maximum distance to\nbitee for each biter",xlab="",ylab="",freq=F,axes=F)
# axis(1,cex.axis=0.7,padj=-0.5)
# axis(2,cex.axis=0.7,padj=0.5)
# mtext("Distance from biter to bitee",1,2,cex=0.7)
# mtext("Density",2,2,cex=0.7)
# box(bty="l")
# lines(seq(0,16500,100),dgamma(seq(0,16500,100),shape=coef(max_dist_gamma)["shape"],rate=coef(max_dist_gamma)["rate"]),
#       col="navy",lwd=1.5)
# lines(rep(qgamma(.975,shape=coef(max_dist_gamma)["shape"],rate=coef(max_dist_gamma)["rate"]), 100), 
#       seq(0,1,length=100), col="navy",lty=1)
# lines(seq(0,16500,100),dlnorm(seq(0,16500,100),meanlog=coef(max_dist_lnorm)["meanlog"],sdlog=coef(max_dist_lnorm)["sdlog"]),
#       col="orange",lwd=1.5,lty=2)
# lines(rep(qlnorm(.975,meanlog=coef(max_dist_lnorm)["meanlog"],sdlog=coef(max_dist_lnorm)["sdlog"]), 100), 
#       seq(0,1,length=100), col="orange",lty=2)
# lines(seq(0,16500,100),dweibull(seq(0,16500,100),shape=coef(max_dist_weibull)["shape"],scale=coef(max_dist_weibull)["scale"]),
#       col="green3",lwd=1.5,lty=3)
# lines(rep(qweibull(.975,shape=coef(max_dist_weibull)["shape"],scale=coef(max_dist_weibull)["scale"]), 100), 
#       seq(0,1,length=100), col="green3",lty=3)
# lines(rep(quantile(rabid_carnivores$max_dist_contact[idx],0.975,na.rm=T), 100), 
#       seq(0,1,length=100), col="lightgrey",lty=4)
# legend("topright",
#        paste(c("gamma","lognormal","weibull"),", AIC=",round(c(max_dist_gamma$aic,max_dist_lnorm$aic,max_dist_weibull$aic)),sep=""),
#        lty=c(1:3),col=c("navy","orange","green3"),bty="n",lwd=1.5,cex=0.7)
# text(15000,0.001,paste("mean=", round(mean(rabid_carnivores$max_dist_contact,na.rm=T),2),sep=""),cex=0.7)
# 
# #  Kernel fitted to all distances
# hist(dist_contact_all$distance, breaks=seq(-0.5,35000.5,500),cex.main=0.9,
#      col="lightgrey",main="Using all distances to\nbitees for each biter",xlab="",ylab="",freq=F,axes=F)
# axis(1,cex.axis=0.7,padj=-0.5)
# axis(2,cex.axis=0.7,padj=0.5)
# mtext("Distance from biter to bitee",1,2,cex=0.7)
# mtext("Density",2,2,cex=0.7)
# box(bty="l")
# lines(seq(0,16500,100),dgamma(seq(0,16500,100),shape=coef(all_dist_gamma)["shape"],rate=coef(all_dist_gamma)["rate"]),
#       col="navy",lwd=1.5)
# lines(rep(qgamma(.975,shape=coef(all_dist_gamma)["shape"],rate=coef(all_dist_gamma)["rate"]), 100), 
#       seq(0,1,length=100), col="navy",lty=1)
# lines(seq(0,16500,100),dlnorm(seq(0,16500,100),meanlog=coef(all_dist_lnorm)["meanlog"],sdlog=coef(all_dist_lnorm)["sdlog"]),
#       col="orange",lwd=1.5,lty=2)
# lines(rep(qlnorm(.975,meanlog=coef(all_dist_lnorm)["meanlog"],sdlog=coef(all_dist_lnorm)["sdlog"]), 100), 
#       seq(0,1,length=100), col="orange",lty=2)
# lines(seq(0,16500,100),dweibull(seq(0,16500,100),shape=coef(all_dist_weibull)["shape"],scale=coef(all_dist_weibull)["scale"]),
#       col="green3",lwd=1.5,lty=3)
# lines(rep(qweibull(.975,shape=coef(all_dist_weibull)["shape"],scale=coef(all_dist_weibull)["scale"]), 100), 
#       seq(0,1,length=100), col="green3",lty=3)
# lines(rep(quantile(dist_contact_all$distance,0.975,na.rm=T), 100), 
#       seq(0,1,length=100), col="lightgrey",lty=4)
# legend("topright",
#        paste(c("gamma","lognormal","weibull"),", AIC=",round(c(all_dist_gamma$aic,all_dist_lnorm$aic,all_dist_weibull$aic)),sep=""),
#        lty=c(1:3),col=c("navy","orange","green3"),bty="n",lwd=1.5,cex=0.7)
# text(15000,0.001,paste("mean=", round(mean(dist_contact_all$distance,na.rm=T),2),sep=""),cex=0.7)
# 
# #  Kernel sum of two max distances
# hist(sim_max_dists_2, breaks=seq(-0.5,158000.5,500),cex.main=0.9,xlim=c(0,29000),
#      col="lightgrey",main="Sum of two maximum distances from\nbiter to bitee, drawn from gamma",xlab="",ylab="",freq=F,axes=F)
# axis(1,cex.axis=0.7,padj=-0.5)
# axis(2,cex.axis=0.7,padj=0.5)
# mtext("Sum of two distances from biter to bitee",1,2,cex=0.7)
# mtext("Density",2,2,cex=0.7)
# box(bty="l")
# lines(seq(0,16500,100),dgamma(seq(0,16500,100),shape=coef(max_dist_gamma_2)["shape"],rate=coef(max_dist_gamma_2)["rate"]),
#       col="navy",lwd=1.5)
# lines(rep(qgamma(.975,shape=coef(max_dist_gamma_2)["shape"],rate=coef(max_dist_gamma_2)["rate"]), 100), 
#       seq(0,1,length=100), col="navy",lty=1)
# lines(rep(quantile(sim_max_dists_2,0.975,na.rm=T), 100), 
#       seq(0,1,length=100), col="lightgrey",lty=4)
# lines(seq(0,16500,100),dlnorm(seq(0,16500,100),meanlog=coef(max_dist_lnorm_2)["meanlog"],sdlog=coef(max_dist_lnorm_2)["sdlog"]),
#       col="orange",lwd=1.5,lty=2)
# lines(rep(qlnorm(.975,meanlog=coef(max_dist_lnorm_2)["meanlog"],sdlog=coef(max_dist_lnorm_2)["sdlog"]), 100), 
#       seq(0,1,length=100), col="orange",lty=2)
# lines(seq(0,16500,100),dweibull(seq(0,16500,100),shape=coef(max_dist_weibull_2)["shape"],scale=coef(max_dist_weibull_2)["scale"]),
#       col="green3",lwd=1.5,lty=3)
# lines(rep(qweibull(.975,shape=coef(max_dist_weibull_2)["shape"],scale=coef(max_dist_weibull_2)["scale"]), 100), 
#       seq(0,1,length=100), col="green3",lty=3)
# lines(rep(quantile(sim_max_dists_2,0.975,na.rm=T), 100), 
#       seq(0,1,length=100), col="lightgrey",lty=4)
# legend("topright",
#        paste(c("gamma","lognormal","weibull"),", AIC=",round(c(max_dist_gamma_2$aic,max_dist_lnorm_2$aic,max_dist_weibull_2$aic)),sep=""),
#        lty=c(1:3),col=c("navy","orange","green3"),bty="n",lwd=1.5,cex=0.7)
# 
# text(20000,0.0002,paste("mean=", round(mean(sim_max_dists_2,na.rm=T),2),sep=""),cex=0.7)
# 
# 
# #  Kernel sum of two distances
# hist(sim_all_dists_2, breaks=seq(-0.5,72000.5,500),cex.main=0.9,xlim=c(0,29000),
#      col="lightgrey",main="Sum of two distances from\nbiter to bitee, drawn from weibull",xlab="",ylab="",freq=F,axes=F)
# axis(1,cex.axis=0.7,padj=-0.5)
# axis(2,cex.axis=0.7,padj=0.5)
# mtext("Sum of two distances from biter to bitee",1,2,cex=0.7)
# mtext("Density",2,2,cex=0.7)
# box(bty="l")
# lines(seq(0,16500,100),dgamma(seq(0,16500,100),shape=coef(all_dist_gamma_2)["shape"],rate=coef(all_dist_gamma_2)["rate"]),
#       col="navy",lwd=1.5)
# lines(rep(qgamma(.975,shape=coef(all_dist_gamma_2)["shape"],rate=coef(all_dist_gamma_2)["rate"]), 100), 
#       seq(0,1,length=100), col="navy",lty=1)
# lines(seq(0,16500,100),dlnorm(seq(0,16500,100),meanlog=coef(all_dist_lnorm_2)["meanlog"],sdlog=coef(all_dist_lnorm_2)["sdlog"]),
#       col="orange",lwd=1.5,lty=2)
# lines(rep(qlnorm(.975,meanlog=coef(all_dist_lnorm_2)["meanlog"],sdlog=coef(all_dist_lnorm_2)["sdlog"]), 100), 
#       seq(0,1,length=100), col="orange",lty=2)
# lines(seq(0,16500,100),dweibull(seq(0,16500,100),shape=coef(all_dist_weibull_2)["shape"],scale=coef(all_dist_weibull_2)["scale"]),
#       col="green3",lwd=1.5,lty=3)
# lines(rep(qweibull(.975,shape=coef(all_dist_weibull_2)["shape"],scale=coef(all_dist_weibull_2)["scale"]), 100), 
#       seq(0,1,length=100), col="green3",lty=3)
# lines(rep(quantile(sim_all_dists_2,0.975,na.rm=T), 100), 
#       seq(0,1,length=100), col="lightgrey",lty=4)
# legend("topright",
#        paste(c("gamma","lognormal","weibull"),", AIC=",round(c(all_dist_gamma_2$aic,all_dist_lnorm_2$aic,all_dist_weibull_2$aic)),sep=""),
#        lty=c(1:3),col=c("navy","orange","green3"),bty="n",lwd=1.5,cex=0.7)
# text(20000,0.0002,paste("mean=", round(mean(sim_all_dists_2,na.rm=T),2),sep=""),cex=0.7)
# 
# 
# #  Kernel sum of two max distances accounting for change of direction
# hist(sim_max_dists_2_direction, breaks=seq(-0.5,130000.5,500),cex.main=0.9,xlim=c(0,26500),
#      col="lightgrey",main="Sum of two max distances from biter to bitee,\naccounting for random direction",xlab="",ylab="",freq=F,axes=F)
# axis(1,cex.axis=0.7,padj=-0.5)
# axis(2,cex.axis=0.7,padj=0.5)
# mtext("Sum of two distances from biter to bitee",1,2,cex=0.7)
# mtext("Density",2,2,cex=0.7)
# box(bty="l")
# lines(seq(0,16500,100),dgamma(seq(0,16500,100),shape=coef(max_dist_gamma_2_direction)["shape"],rate=coef(max_dist_gamma_2_direction)["rate"]),
#       col="navy",lwd=1.5)
# lines(rep(qgamma(.975,shape=coef(max_dist_gamma_2_direction)["shape"],rate=coef(max_dist_gamma_2_direction)["rate"]), 100), 
#       seq(0,1,length=100), col="navy",lty=1)
# lines(seq(0,16500,100),dlnorm(seq(0,16500,100),meanlog=coef(max_dist_lnorm_2_direction)["meanlog"],sdlog=coef(max_dist_lnorm_2_direction)["sdlog"]),
#       col="orange",lwd=1.5,lty=2)
# lines(rep(qlnorm(.975,meanlog=coef(max_dist_lnorm_2_direction)["meanlog"],sdlog=coef(max_dist_lnorm_2_direction)["sdlog"]), 100), 
#       seq(0,1,length=100), col="orange",lty=2)
# lines(seq(0,16500,100),dweibull(seq(0,16500,100),shape=coef(max_dist_weibull_2_direction)["shape"],scale=coef(max_dist_weibull_2_direction)["scale"]),
#       col="green3",lwd=1.5,lty=3)
# lines(rep(qweibull(.975,shape=coef(max_dist_weibull_2_direction)["shape"],scale=coef(max_dist_weibull_2_direction)["scale"]), 100), 
#       seq(0,1,length=100), col="green3",lty=3)
# lines(rep(quantile(sim_max_dists_2_direction,0.975,na.rm=T), 100), 
#       seq(0,1,length=100), col="lightgrey",lty=4)
# legend("topright",
#        paste(c("gamma","lognormal","weibull"),", AIC=",round(c(max_dist_gamma_2_direction$aic,max_dist_lnorm_2_direction$aic,max_dist_weibull_2_direction$aic)),sep=""),
#        lty=c(1:3),col=c("navy","orange","green3"),bty="n",lwd=1.5,cex=0.7)
# text(20000,0.00025,paste("mean=", round(mean(sim_max_dists_2_direction,na.rm=T),2),sep=""),cex=0.7)
# 
# 
# #  Kernel sum of two distances accounting for change of direction
# hist(sim_dists_2_direction, breaks=seq(-0.5,66000.5,500),cex.main=0.9,xlim=c(0,26500),
#      col="lightgrey",main="Sum of two distances from biter to bitee,\naccounting for random direction",xlab="",ylab="",freq=F,axes=F)
# axis(1,cex.axis=0.7,padj=-0.5)
# axis(2,cex.axis=0.7,padj=0.5)
# mtext("Sum of two distances from biter to bitee",1,2,cex=0.7)
# mtext("Density",2,2,cex=0.7)
# box(bty="l")
# lines(seq(0,16500,100),dgamma(seq(0,16500,100),shape=coef(all_dist_gamma_2_direction)["shape"],rate=coef(all_dist_gamma_2_direction)["rate"]),
#       col="navy",lwd=1.5)
# lines(rep(qgamma(.975,shape=coef(all_dist_gamma_2_direction)["shape"],rate=coef(all_dist_gamma_2_direction)["rate"]), 100), 
#       seq(0,1,length=100), col="navy",lty=1)
# lines(seq(0,16500,100),dlnorm(seq(0,16500,100),meanlog=coef(all_dist_lnorm_2_direction)["meanlog"],sdlog=coef(all_dist_lnorm_2_direction)["sdlog"]),
#       col="orange",lwd=1.5,lty=2)
# lines(rep(qlnorm(.975,meanlog=coef(all_dist_lnorm_2_direction)["meanlog"],sdlog=coef(all_dist_lnorm_2_direction)["sdlog"]), 100), 
#       seq(0,1,length=100), col="orange",lty=2)
# lines(seq(0,16500,100),dweibull(seq(0,16500,100),shape=coef(all_dist_weibull_2_direction)["shape"],scale=coef(all_dist_weibull_2_direction)["scale"]),
#       col="green3",lwd=1.5,lty=3)
# lines(rep(qweibull(.975,shape=coef(all_dist_weibull_2_direction)["shape"],scale=coef(all_dist_weibull_2_direction)["scale"]), 100), 
#       seq(0,1,length=100), col="green3",lty=3)
# lines(rep(quantile(sim_dists_2_direction,0.975,na.rm=T), 100), 
#       seq(0,1,length=100), col="lightgrey",lty=4)
# legend("topright",
#        paste(c("gamma","lognormal","weibull"),", AIC=",round(c(all_dist_gamma_2_direction$aic,all_dist_lnorm_2_direction$aic,all_dist_weibull_2_direction$aic)),sep=""),
#        lty=c(1:3),col=c("navy","orange","green3"),bty="n",lwd=1.5,cex=0.7)
# text(20000,0.00025,paste("mean=", round(mean(sim_dists_2_direction,na.rm=T),2),sep=""),cex=0.7)


# dev.off()

