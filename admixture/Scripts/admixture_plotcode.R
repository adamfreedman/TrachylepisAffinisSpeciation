popnames=c("Lope-gallery","Lope-forest","Kribi","Mbongwana","Palm d'Or","Malimba","Mapenja","Konye","Mundemba","Nguti","Doume","Mbakaou","Bounou","Bazzama")

#
#K7 BEST
data<-read.table("withids_oneperrad_trachylepis.7.Q",header=FALSE)
colnames(data)<-c("id","pop","c0","c1","c2","c3","c4","c5","c6")
popcounts<-c()
for(i in seq(1,14)){popcounts<-c(popcounts,length(data$pop[data$pop==i]))}
sumcounts<-c()
for(i in seq(2,14)){sumcounts<-c(sumcounts,sum(popcounts)-sum(popcounts[i:14]))}
labelpos<-c(c(0),sumcounts,c(208))
meanpos<-c()
for(i in seq(1,14)){meanpos<-c(meanpos,(labelpos[i]+(labelpos[i+1]-labelpos[i])/2))}

sortdata<-data[order(data$pop),]
pdf(file="Trachylepis_admixtureplotK7_nofilt_2020.12.08.pdf",width=7,height=4)
barplot(t(as.matrix(sortdata[,3:9])),col=c("darkorchid","magenta","red","cyan","dodgerblue","gold","blue"),xaxt="n",ylab="Ancestry",border=NA,space=0)
for(i in seq(1,13)){segments(sumcounts[i]+0,0,sumcounts[i]+0,1,col="black",lwd=1.5)}
mtext(popnames,side=1,at=meanpos,las=2,line=0)
dev.off()


#K2
data<-read.table("withids_oneperrad_trachylepis.2.Q",header=FALSE)
colnames(data)<-c("id","pop","c0","c1")
####
#popcounts<-c()
#for(i in seq(1,14)){popcounts<-c(popcounts,length(data$pop[data$pop==i]))}
#sumcounts<-c()
#for(i in seq(2,14)){sumcounts<-c(sumcounts,sum(popcounts)-sum(popcounts[i:14]))}
#meanpos<-c()
#for(i in seq(1,14)){meanpos<-c(meanpos,(labelpos[i]+(labelpos[i+1]-labelpos[i])/2))}
####
colnames(data)<-c("id","pop","c0","c1")
sortdata<-data[order(data$pop),]
pdf(file="Trachylepis_admixtureplotK2_nofilt_2020.12.08.pdf",width=7,height=4)
barplot(t(as.matrix(sortdata[,3:4])),col=rainbow(2),xaxt="n",ylab="Ancestry",border=NA,space=0)
for(i in seq(1,13)){segments(sumcounts[i]+0,0,sumcounts[i]+0,1,col="black",lwd=1.5)}
mtext(popnames,side=1,at=meanpos,las=2,line=0)
dev.off()



#K3

data<-read.table("withids_oneperrad_trachylepis.3.Q",header=FALSE)
colnames(data)<-c("id","pop","c0","c1","c2")
sortdata<-data[order(data$pop),]

pdf(file="Trachylepis_admixtureplotK3_nofilt_2020.12.08.pdf",width=7,height=4)
barplot(t(as.matrix(sortdata[,3:5])),col=c("#00FFFFFF","#FF0000FF","#0000FFFF"),xaxt="n",ylab="Ancestry",border=NA,space=0)
for(i in seq(1,13)){segments(sumcounts[i],0,sumcounts[i],1,col="black",lwd=1.5)}

mtext(popnames,side=1,at=meanpos,las=2,line=0)
dev.off()



#K4

data<-read.table("withids_oneperrad_trachylepis.4.Q",header=FALSE)
colnames(data)<-c("id","pop","c0","c1","c2","c4")
sortdata<-data[order(data$pop),]

pdf(file="Trachylepis_admixtureplotK4_nofilt_2020.12.08.pdf",width=7,height=4)
barplot(t(as.matrix(sortdata[,3:6])),col=c("blue","gold","cyan","red"),xaxt="n",ylab="Ancestry",border=NA,space=0)
for(i in seq(1,13)){segments(sumcounts[i],0,sumcounts[i],1,col="black",lwd=1.5)}
mtext(popnames,side=1,at=meanpos,las=2,line=0)
dev.off()


#K5

data<-read.table("withids_oneperrad_trachylepis.5.Q",header=FALSE)
colnames(data)<-c("id","pop","c0","c1","c2","c4","c5")

sortdata<-data[order(data$pop),]

pdf(file="Trachylepis_admixtureplotK5_nofilt_2020.12.08.pdf",width=7,height=4)
barplot(t(as.matrix(sortdata[,3:7])),col=c("red","dodgerblue","cyan","blue","gold"),xaxt="n",ylab="Ancestry",border=NA,space=0)
for(i in seq(1,13)){segments(sumcounts[i],0,sumcounts[i],1,col="black",lwd=1.5)}
mtext(popnames,side=1,at=meanpos,las=2,line=0)
dev.off()



#K6

data<-read.table("withids_oneperrad_trachylepis.6.Q",header=FALSE)
colnames(data)<-c("id","pop","c0","c1","c2","c3","c4","c5")

sortdata<-data[order(data$pop),]

pdf(file="Trachylepis_admixtureplotK6_nofilt_2020.12.08.pdf",width=7,height=4)
barplot(t(as.matrix(sortdata[,3:8])),col=c("red","dodgerblue","blue","magenta","gold","cyan"),xaxt="n",ylab="Ancestry",border=NA,space=0)
for(i in seq(1,13)){segments(sumcounts[i]+0.5,0,sumcounts[i]+0.5,1,col="black",lwd=1)}
mtext(popnames,side=1,at=meanpos,las=2,line=0)
dev.off()


