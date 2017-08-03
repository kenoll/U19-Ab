setwd("~/Dropbox/Heise/U19-Ab/weight/chr1_microsatellite")
prect=read.csv("microsatellite_queries/preCTCTT_kmer_query1_1_singles.csv")
row.names(prect)<-prect$strain
prect=prect[-(1:2)]

prect[prect > 1] = 5

k <- kmeans(prect, 7)
df=prect
dfc <- cbind(df, cluster=k$cluster)

write.csv(dfc,"microsatellite_queries/preCTCTT_kmer_query1_1_clusters.csv")

dfc=read.csv("microsatellite_queries/preCTCTT_kmer_query1_1_clusters.csv")
colnames(dfc)[1]="strain"
score=read.csv("~/Dropbox/Heise/U19-Ab/weight/weight_loss/qtls/chr1_38-41_norecomb_haplotypes.csv")

pheno=dfc
pheno=merge(pheno,score)
write.csv(pheno,"preCTCTT_cluster_scores.csv")

ggplot(pheno,aes(x=founders,y=cluster))+
  geom_point(aes(color=status))+
  geom_text(aes(label=strain),position="jitter",size=2)

?ggplot

dfc$idsort <- dfc$id[order(dfc$cluster)]
dfc$idsort <- order(dfc$idsort)

dfm <- melt(dfc, id.vars=c("id", "idsort"))

ggplot(dfm, aes(x=variable, y=idsort)) + geom_tile(aes(fill=value))


ggplot(prectm, aes(variable, founders)) + geom_tile() +
  scale_fill_gradient(low = "white", high = "steelblue")



##mat
prect=as.matrix(prect)
prect=t(prect)

heatmap(prect, Colv=F, scale='none')

