## sims ##
no_exp<-read.table("ForVsEcoBestModNoExpansion.weir.fst",header=TRUE)
bestmod<-read.table("ForVsEcoBestModel.weir.fst",header=TRUE)
no_divbottle<-read.table("ForVsEcoBestModelNoDivBottleneck.weir.fst",header=TRUE)
threeepoch<-read.table("anc_asym_mig_size_3epoch_earlysecondarycontact_simfromsfs.fst",header=TRUE)

## data ##
foreco<-read.table("maf05_allforestVsecotone.weir.fst.weir.fst",header=TRUE)
foreco_fromsfs<-read.table("ForVsEcoDataFromSFS.weir.fst",header=TRUE)

### histograms ###
no_exp_hist<-hist(no_exp$WEIR_AND_COCKERHAM_FST,breaks=seq(-.05,1,.05),plot=F)
foreco_emp_hist<-hist(foreco$WEIR_AND_COCKERHAM_FST,breaks=seq(-.05,1,.05),plot=F)
best_hist<-hist(bestmod$WEIR_AND_COCKERHAM_FST,breaks=seq(-.05,1,.05),plot=F)
foreco_fromsfs_hist<-hist(foreco_fromsfs$WEIR_AND_COCKERHAM_FST,breaks=seq(-.05,1,.05),plot=F)
nodivbottle_hist<-hist(no_divbottle$WEIR_AND_COCKERHAM_FST,breaks=seq(-.05,1,.05),plot=F)
three_epoch_hist<-hist(threeepoch$WEIR_AND_COCKERHAM_FST,breaks=seq(-.05,1,.05),plot=F)


pdf(file="GenomicFstDistributions_ForestVsEcotone.pdf",width=4,height=4)
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1),ylim=c(0,4))
lines(no_exp_hist$mids,no_exp_hist$density,col="black",lwd=2,lty=2)
lines(foreco_emp_hist$mids,foreco_emp_hist$density,col="maroon",lwd=2,lty=2)
lines(best_hist$mids,best_hist$density,col="dodgerblue",lwd=2,lty=2)
lines(foreco_fromsfs_hist$mids,foreco_fromsfs_hist$density,col="yellow3",lwd=2,lty=2)
lines(nodivbottle_hist$mids,nodivbottle_hist$density,col="darkorange3",lwd=2,lty=2)
lines(three_epoch_hist$mids,three_epoch_hist$density,col="blue",lwd=2,lty=2)


mtext(expression(paste(F[ST],": forest vs. ecotone")),side=1,line=2.3)
mtext("Density",side=2,line=2.3)
legend("topright",legend=c("Data","Data: from SFS","Best Model: ancient gene flow + expansion","Model: no expansion","Model: no divergence bottleneck", "Model: 3-epoch"),cex=0.6,col=c("maroon","yellow3","dodgerblue","black","darkorange3","blue"),lty=1,lwd=2)
dev.off()
