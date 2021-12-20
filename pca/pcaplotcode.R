pcadata<-read.table("wpoplabels_taff_smartpca_10pcs_Rtable.txt",header=TRUE)
# pc1 vs pc2
pdf("taffinis_snppca_pc1v2_2020.08.14.pdf",width=5,height=5)
plot(1, type="n", xaxt="n",yaxt="n",xlab="", ylab="", xlim=c(min(pcadata$pc1),max(pcadata$pc1)), ylim=c(min(pcadata$pc2),max(pcadata$pc2)))
axis(side=1,at=c(0,-0.05,0,0.05,0.10),cex.axis=0.8,padj=-1)
axis(side=2,at=c(-0.15,-0.10,-0.05,0,0.05),cex.axis=0.8,pad=1)
mtext("PC1 (30.0%)",side=1,line=2,cex=1)
mtext("PC2 (17.8%)",side=2,line=2,cex=1)

# MBK
#points(pcadata$pc1[pcadata$pop==12],pcadata$pc2[pcadata$pop==12],col=rgb(2,146,255,150,max=255),pch=16)
points(pcadata$pc1[pcadata$pop==12],pcadata$pc2[pcadata$pop==12],col=rainbow(7)[5],pch=3)

#DOU
#points(pcadata$pc1[pcadata$pop==11],pcadata$pc2[pcadata$pop==11],col=rgb(2,146,255,150,max=255),pch=1)
points(pcadata$pc1[pcadata$pop==11],pcadata$pc2[pcadata$pop==11],col=rainbow(7)[5],pch=3)

# BOUNOU
#points(pcadata$pc1[pcadata$pop==13],pcadata$pc2[pcadata$pop==13],col=rgb(0,100,0,150,max=255),pch=3)
points(pcadata$pc1[pcadata$pop==13],pcadata$pc2[pcadata$pop==13],col=rainbow(7)[4],pch=3)

# NGUTI
#points(pcadata$pc1[pcadata$pop==10],pcadata$pc2[pcadata$pop==10],col=rgb(73,0,255,150,max=255),pch=16)
points(pcadata$pc1[pcadata$pop==10],pcadata$pc2[pcadata$pop==10],col=rainbow(7)[6],pch=3)

# MUNDEMBA
points(pcadata$pc1[pcadata$pop==9],pcadata$pc2[pcadata$pop==9],col=rainbow(7)[6],pch=3)

#KONYE
points(pcadata$pc1[pcadata$pop==8],pcadata$pc2[pcadata$pop==8],col=rainbow(7)[6],pch=3)

# BAZ
points(pcadata$pc1[pcadata$pop==14],pcadata$pc2[pcadata$pop==14],col=rainbow(7)[4],pch=3)

#KRI
points(pcadata$pc1[pcadata$pop==3],pcadata$pc2[pcadata$pop==3],col=rainbow(7)[7],pch=3)

#MBW
points(pcadata$pc1[pcadata$pop==4],pcadata$pc2[pcadata$pop==4],col=rainbow(7)[7],pch=3)

#LOP-GAL
points(pcadata$pc1[pcadata$pop==1],pcadata$pc2[pcadata$pop==1],col=rainbow(7)[3],pch=3)

#LOP-FOR
points(pcadata$pc1[pcadata$pop==2],pcadata$pc2[pcadata$pop==2],col=rainbow(7)[3],pch=3)

#MPJ
points(pcadata$pc1[pcadata$pop==7],pcadata$pc2[pcadata$pop==7],col=rainbow(7)[2],pch=3)

#POR
points(pcadata$pc1[pcadata$pop==5],pcadata$pc2[pcadata$pop==5],col=rgb(60,60,60,100,max=255),pch=3)

#MAL

points(pcadata$pc1[pcadata$pop==6],pcadata$pc2[pcadata$pop==6],col=rainbow(7)[1],pch=3)

legend("bottomright",legend=c("Lope NP: Gabon","South of Sanaga","South of Sanaga: admixed","North of Sanaga","Montane forest","Southwest provice","Western Ecotone","Eastern ecotone"),pch=3,col=c(rainbow(7)[3],rainbow(7)[7],rgb(60,60,60,100,max=255),rainbow(7)[1],rainbow(7)[2],rainbow(7)[6],rainbow(7)[5],rainbow(7)[4]),cex=0.8)
dev.off()


