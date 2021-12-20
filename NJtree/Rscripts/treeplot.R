ids<-read.table("wcols_ids_renumbered.txt",header=TRUE)
ibsdist<-read.table("plink.mdist",header=FALSE)
ibsdist<-as.matrix(ibsdist)
rownames(ibsdist)<-ids$id
colnames(ibsdist)<-ids$id
tr<-nj(ibsdist)
pdf("NJtree_withlabels.pdf")
plot.phylo(tr,"unrooted",tip.color=ids$color,lab4ut="axial")
dev.off()

pdf("NJtree_withoutlabels.pdf")
plot.phylo(tr,"unrooted",show.tip.label=FALSE,lab4ut="axial")
dev.off()
