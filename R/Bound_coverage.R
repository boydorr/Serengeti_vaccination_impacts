
pdf("Figs/Bound_coverage_function.pdf",width=4.5, height=3.4)

par(mfrow=c(1,1),mar=c(5,5,1,1))

C<-seq(0,1.5,0.01)
plot(c(0,max(C)),c(0,max(C)),lty=2,type="l",bty="l",ylim=c(0,1.2),xlim=c(0,max(C)),ylab="b(x)",xlab="x",cex.axis=0.8)
lines(c(0,max(C)),c(1,1),lty=2)
lines(c(1,1),c(0,max(C)),lty=2)
cols <- brewer.pal(5,"Set2")
# lines(c(0,max(C)),c(0.1,0.1),lty=2)
a<-1; lines(C,C/(1+C^a)^(1/a),col=cols[1])
a<-3; lines(C,C/(1+C^a)^(1/a),col=cols[2])
a<-4; lines(C,C/(1+C^a)^(1/a),col=cols[3])
a<-5; lines(C,C/(1+C^a)^(1/a),col=cols[4])
a<-9; lines(C,C/(1+C^a)^(1/a),col=cols[5])
a<-6; lines(C,C/(1+C^a)^(1/a),col=1,lwd=2)
legend(0.03,1.225,paste("a =",c(1,3,4,5,6,9)),col=c(cols[1:4],1,cols[5]),lwd=c(1,1,1,1,2,1),cex=0.8,bg="white")

dev.off()

a<-6
0.3/(1+0.3^a)^(1/a)
0.5/(1+0.5^a)^(1/a)
0.6/(1+0.6^a)^(1/a)
0.7/(1+0.7^a)^(1/a)
1/(1+1^a)^(1/a)
2/(1+2^a)^(1/a)
3/(1+3^a)^(1/a)



