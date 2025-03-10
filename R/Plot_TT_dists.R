rm(list=ls())

set.seed(0)

# Read in serial intervals and serial interval distribution
SI_params <- read.csv("output/SI_params.csv")
SIs <- read.table("output/serial_intervals.csv")[,1]
length(SIs)

# Read in distances and distance kernel
DK_params <- read.csv("output/DK_params.csv")
distances <- read.table("output/distances_between_case_contacts.csv")[,1]
length(distances)

# Get probability for different serial intervals
SI_seq <- seq(0,max(SIs),1)
fitted_values_SI <- dlnorm(SI_seq,meanlog = SI_params$SI_ml,sdlog = SI_params$SI_sdlog)

# SI mean and CIs
SI_means <- rep(NA,10000)
for(i in 1:10000){
  SI_means[i] <- mean(sample(SIs,size=length(SIs),replace = T))}
SI_CI <- quantile(SI_means,c(0.025,0.975))

# Get probability for different distances
distance_seq <- seq(0,max(distances),5)
row <- which(DK_params$kernel_type=="DK" & DK_params$dist=="weibull")
DK_shape = DK_params$par1est[row]; DK_scale = DK_params$par2est[row] 
fitted_values_distance <- dweibull(distance_seq,shape=DK_shape,scale=DK_scale)

# Distance mean and CIs
distance_means <- rep(NA,10000)
for(i in 1:10000){
  distance_means[i] <- mean(sample(distances,size=length(distances),replace = T))}
distance_CI <- quantile(distance_means,c(0.025,0.975))


#Plot
cex.axis=0.75
cex.lab=0.9
cex.main=0.9
pdf("Figs/DK_SI_dist.pdf",width=7, height=3.5)
par(mfrow=c(1,2))
par(mar=c(3,2.3,3,0.9))
interval <- 4
hist(SIs,breaks=seq(-0.5,max(SIs)+interval,interval),
     ylim=c(0,max(fitted_values_SI)),
     xlim=c(0,max(SIs)),
     border=F,main="Serial Interval Distribution",freq=F,xlab="",ylab="",axes=F,cex.main=cex.main)
axis(2,cex.axis=cex.axis,padj=1)
axis(1,cex.axis=cex.axis,padj=-1.5)
box(bty="l")
mtext("Density",side=2,line=1.5,cex=cex.lab)
mtext("Serial Interval (days)",side=1,line=2,cex=cex.lab)
box(bty="l")
lines(fitted_values_SI~SI_seq,lwd=2,col="navy")
legend(-0.025*max(SIs),0.85*max(fitted_values_SI[-1]),
       c(paste("data (mean=", round(mean(SIs),1),"(95%CI: ", paste0(round(SI_CI,1),collapse = "-"),"), n=", length(SIs),"),",sep=""),
         paste("lognormal fit (meanlog=", round(SI_params$SI_ml,2),", sdlog=",round(SI_params$SI_sdlog,2),")",sep="")),
       text.col = c("grey40","navy"),bty="n",y.intersp = 1.5,cex=cex.axis)
legend("topleft","A",text.font=2,bty="n",cex=1.1)

interval <-200
hist(distances,breaks=seq(-0.5,max(distances)+interval,interval),
     ylim=c(0,max(fitted_values_distance[-1])),
     xlim=c(0,max(distances)),
     border=F,main="Distance Kernel",freq=F,xlab="",ylab="",axes=F,cex.main=cex.main)
axis(2,cex.axis=cex.axis,padj=1)
axis(1,cex.axis=cex.axis,padj=-1.5)
box(bty="l")
mtext("Density",side=2,line=1.5,cex=cex.lab)
mtext("Distance (km)",side=1,line=2,cex=cex.lab)
box(bty="l")
lines(fitted_values_distance~distance_seq,lwd=2,col="navy")
legend(-0.025*max(distances),0.85*max(fitted_values_distance[-1]),
       c(paste("data (mean=", round(mean(distances)),"(95%CI: ", paste0(round(distance_CI),collapse = "-"),"), n=", length(distances),"),",sep=""),
         paste("weibull fit (shape=", round(DK_shape,2),", scale=",round(DK_scale,2),")",sep="")),
       text.col = c("grey40","navy"),bty="n",y.intersp = 1.5,cex=cex.axis)
legend("topleft","B",text.font=2,bty="n",cex=1.1)

dev.off()



