eems1<-read.table("oneperrad_trachylepis-EEMS-nDemes200-1/mcmcpilogl.txt",header=FALSE)
colnames(eems1)<-c("prior","likelihood")
eems1$logposterior<-eems1$prior+eems1$likelihood

eems2<-read.table("oneperrad_trachylepis-EEMS-nDemes200-2/mcmcpilogl.txt",header=FALSE)
colnames(eems2)<-c("prior","likelihood")
eems2$logposterior<-eems2$prior+eems2$likelihood

eems3<-read.table("oneperrad_trachylepis-EEMS-nDemes200-3/mcmcpilogl.txt",header=FALSE)
colnames(eems3)<-c("prior","likelihood")
eems3$logposterior<-eems3$prior+eems3$likelihood

eems4<-read.table("oneperrad_trachylepis-EEMS-nDemes200-4/mcmcpilogl.txt",header=FALSE)
colnames(eems4)<-c("prior","likelihood")
eems4$logposterior<-eems4$prior+eems4$likelihood

eems5<-read.table("oneperrad_trachylepis-EEMS-nDemes200-5/mcmcpilogl.txt",header=FALSE)
colnames(eems5)<-c("prior","likelihood")
eems5$logposterior<-eems5$prior+eems5$likelihood

pdf("EeemsChainMixing.pdf",width=5,height=3)
plot(1, type="n", xlab="MCMC sample", ylab="log posterior", xlim=c(0,100), ylim=c(80100,80300))
lines(seq(1,100),eems1$logposterior,col="firebrick",lwd=2)
lines(seq(1,100),eems2$logposterior,col="blue",lwd=2)
lines(seq(1,100),eems3$logposterior,col="darkorchid",lwd=2)
lines(seq(1,100),eems4$logposterior,col="goldenrod",lwd=2)
lines(seq(1,100),eems5$logposterior,col="black",lwd=2)
dev.off()
