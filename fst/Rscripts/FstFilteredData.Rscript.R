library(lattice)
fst<-read.csv("taffinis_pairwiseFstdata_filtered_withdiagonal.csv",header=TRUE)
an <- with(fst, sort(unique(c(as.character(pop1),as.character(pop2)))))
M <- array(0, c(length(an), length(an)), list(an, an))
i <- match(fst$pop1, an)
j <- match(fst$pop2, an)
M[cbind(i,j)] <- M[cbind(j,i)] <- fst$fst
slatkinM<-M/(1-M) 
col.order <- c("lopegallery","lopeforest","kribi","mbongwana","palmdor","malimba","mapenja","konye","mundemba","nguti","doume","mbakaou","bounou","bazzama")
sort1<-M[,col.order]
sort2<-sort1[col.order,]
new.palette=colorRampPalette(c("navy", "white", "firebrick3"))(100)
new.labels<-c("LP1","LP2","KRI","MBW","POR","MAL","MPJ","KON","MUN","NGU","DOU","MBK","BOU","BAZ")

# Weir and Cockerham Fst matrix plot #
pdf(file="taffinis_Fst_maf05.pdf",width=4.5,height=4.5)
levelplot(sort2,col.regions=new.palette,xlab="",ylab="",scales = list(labels=new.labels,x=list(rot=90),tck = c(1,0)))
dev.off()

# Slatkin's Linearlized Fst matrix plot #
sort3<-slatkinM[,col.order]
sort4<-sort3[col.order,]
pdf(file="taffinis_slatkinFst_maf05.pdf",width=4.5,height=4.5)
levelplot(sort4,col.regions=new.palette,xlab="",ylab="",scales = list(labels=new.labels,x=list(rot=90),tck = c(1,0)))
dev.off()