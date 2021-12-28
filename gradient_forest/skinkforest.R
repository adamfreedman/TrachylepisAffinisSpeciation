library(gradientForest)
library(data.table)
library(raster)
library(gstat)

setwd("~/Desktop/R/skinks")
Eskink=read.table("Eskink.txt",header=T,sep="") #Environmental data at sites
Gskink=read.table("Gskink2.txt",header=F,sep="") #Genomic data at sites
skinksites=Eskink[,1:2] #geographic coordinates at sites


Eskink2=Eskink[,3:18] #setting environmental predictors
preds <- colnames(Eskink2) #shorthand for predcitors
specs <- colnames(Gskink) #shorthand for response

nSites <- dim(Gskink)[1]
nSpecs <- dim(Gskink)[2]
# set depth of conditional permutation
lev <- floor(log2(nSites*0.368/2))
lev

#the forest is commented out and loaded from a save here, it takes a bit of time to run
#skinkforest=gradientForest(cbind(Eskink2,Gskink), predictor.vars=preds, response.vars=specs, ntree=600, transform = NULL, compact=F,nbin=100, maxLevel=1,trace=T)
skinkforest=readRDS(file = "skinkforest2.RData")


#randomizations of predictors and response pairing to detect signal apart from random associations

realtotal=skinkforest$species.pos.rsq #observed number of SNPs with positive R squared
realaverage=sum(skinkforest$result)/realtotal #observed average R squared for those SNPs

# here we run 100 randomizations using the observed predictors and responses, but randomly paired
for (i in 1:100){
  
  predsR=Eskink2[sample(nrow(Eskink2)),]
  
  skinkforestR=gradientForest(cbind(predsR,Gskink), predictor.vars=colnames(predsR), response.vars=colnames(Gskink), ntree=200, transform = NULL, compact=T,nbin=100, maxLevel=1,trace=F)
  
  randtotal=skinkforestR$species.pos.rsq
  randaverage=sum(skinkforestR$result)/randtotal
  
  write.table(randtotal,file=paste("randtotal",sep=""),row.names=FALSE,col.names=FALSE,append=T)
  write.table(randaverage,file=paste("randaverage",sep=""),row.names=FALSE,col.names=FALSE,append=T)
}

randa=read.csv("randaverage",sep="",header=T)
randt=read.csv("randtotal",sep="",header=T)

#and here we draw those randomizations and place red lines where we observed our actual forest results (just copied from our real output above)

hist(randa$randa,xlim=c(0,0.5),main="Average R-squared of Random Gradient Forests",xlab="Average R-squared of SNPs")
abline(v=0.381531,col="red")
hist(randt$randt,xlim=c(0,12000),main="Total Correlated SNPs of Random Gradient Forests",xlab="Total SNPs of R-squared > 1")
abline(v=9765,col="red")

#using our forest to predict across the range of skinks
skinkgrid=read.table("skinksm100k.csv",header=T,sep=",") #random 100K points chosen with top predictor values included

#predict overall importance
plot.gradientForest(skinkforest,plot.type="O")

allpredictors=names(importance(skinkforest))

#split density plots
plot(skinkforest, plot.type="S", imp.vars=allpredictors, leg.posn="topright", cex.legend=0.4, cex.axis=0.6, cex.lab=0.7, line.ylab=0.9, par.args=list(mgp=c(1.5, 0.5, 0), mar=c(3.1,1.5,0.1,1)))

#species cumulative plot
plot.gradientForest(skinkforest, plot.type="Cumulative.Importance", imp.vars=allpredictors, show.overall=T, legend=T,common.scale=T,leg.posn="topleft", leg.nspecies=5, cex.lab=0.7, cex.legend=0.4, cex.axis=0.6, line.ylab=0.9, par.args=list(mgp=c(1.5, 0.5, 0), mar=c(3.1,1.5,0.1,1),omi=c(0,0.3,0,0)))

#predictor cumulatives
plot(skinkforest, plot.type="C", imp.vars=allpredictors, show.species=F, common.scale=T, cex.axis=0.6, cex.lab=0.7, line.ylab=0.9, par.args=list(mgp=c(1.5, 0.5, 0), mar=c(2.5,1.0,0.1,0.5), omi=c(0,0.3,0,0)))

#R squared plot
plot(skinkforest, plot.type="P", show.names=F, horizontal=F, cex.axis=1, cex.labels=0.7, line=2.5)

#transform grid and environmental predictors
predictors <- names(importance(skinkforest)[1:3])
tskinkgrid=cbind(skinkgrid[,c("Long","Lat")], predict(skinkforest,skinkgrid[,predictors]))
Trns_site <- predict(skinkforest)

#conduct PC analyses on top three predictors
PCs=prcomp(tskinkgrid[,3:5])
sgn <- sign(PCs$rotation["quickscatstd",])
PCs$rotation <- sweep(PCs$rotation,2,sgn,"*")
PCs$x <- sweep(PCs$x,2,sgn,"*")
# set up a colour palette for the mapping
a1 <- PCs$x[,1]
a2 <- PCs$x[,2]
a3 <- PCs$x[,3]
r <- a1+a2
g <- -a2
b <- a3+a2-a1
r <- (r-min(r)) / (max(r)-min(r)) * 255
g <- (g-min(g)) / (max(g)-min(g)) * 255
b <- (b-min(b)) / (max(b)-min(b)) * 255

nvs <- dim(PCs$rotation)[1] # number of variables
vec <- c("bio2","bio12","quickscatstd") 
lv <- length(vec)
vind <- rownames(PCs$rotation) %in% vec
# choose a scaling factor to plot the vectors over the grid
scal <- 20
xrng <- range(PCs$x[,1], PCs$rotation[,1]/scal)*1.1
yrng <- range(PCs$x[,2], PCs$rotation[,2]/scal)*1.1
plot((PCs$x[,1:2]), xlim=xrng, ylim=yrng, pch=".", cex=7, col=rgb(r,g,b, max = 255), asp=1)
# plot the other predictors with "+"
points(PCs$rotation[! vind,1:2]/scal, pch="+")  
# plot the chosen predictors as arrows
arrows(rep(0,lv), rep(0,lv), PCs$rotation[,1]/scal, PCs$rotation[,2]/scal, length = 0.1)
jit <- 0.0015
text(PCs$rotation[vec,1]/scal+jit*sign(PCs$rotation[vec,1]), PCs$rotation[vec,2]/scal+jit*sign(PCs$rotation[vec,2]), labels = vec)

# first predict the PCs for the transformed site data
PCsites <- predict(PCs,Trns_site[,predictors])
# plot all the sites as points on the biplot
points(PCsites[,1:2])

#plot these in across species range
skink.pred <- predict(skinkforest, skinkgrid[,predictors])
plot(tskinkgrid[,c("Long","Lat")],pch=15,cex=0.5,asp=1,col=rgb(r,g,b, max=255),main="SNP turnover (n=14,991) across the range of T. affinis in Central Africa")
points(Eskink[,1:2])
