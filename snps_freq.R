library(plyr)
library(ggplot2)
library(dplyr)

snps=read.csv("~/Dropbox/Heise/ELISA Antibody/qtls/d10/pre-feb2017/snps/snps_11-69082450-72556480.csv")


snps.count=count(snps$Gene)
colnames(snps.count)=c("Gene","Freq")

barplot(snps.count$Freq, space=0, names = snps.count$Gene, xlab='Gene', ylab = 'Freq')

snps.count[which.max(snps.count$Freq),]
subset(snps.count,snps.count$Freq>500)
median(snps.count$Freq)



colnames(snps)

### genes for which WSB has a missense ###
gene.list= c("Wscd1","Zbtb4","Zfp3","Tnk1","Rpain","Rpain","Slc16a13","Kctd11","Nlrp1b","Nlrp1a")

snps=snps[which(snps$Gene %in% gene.list),]
snps.count=count(snps$Gene)
colnames(snps.count)=c("Gene","Freq")
snps.count







#missense only
snps=read.csv("~/Dropbox/Heise/ELISA Antibody/qtls/d10/pre-feb2017/snps/snps_11-69082450-72556480_missenseonly.csv")

snps=snps[which(snps$Gene %in% gene.list),]
colnames(snps.count)=c("Gene","Freq")
snps.count=count(snps$Gene)
