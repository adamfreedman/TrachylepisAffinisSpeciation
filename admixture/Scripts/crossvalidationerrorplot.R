library(ggplot2)
setwd("/Users/adamfreedman/Dropbox/Trachylepis/admixture")
data<-read.table("CrossValidationErrorTable.tsv",header=TRUE)
p<-ggplot(data=data, aes(x=K, y=CVerror)) + geom_bar(stat="identity") + ylab("Cross-validation error") + scale_x_continuous(breaks=c(2,3,4,5,6,7,8,9,10,11,12))

ggsave("trachylepis_admixture_CVerror.pdf",width=5,height=3)