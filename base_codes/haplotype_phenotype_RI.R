setwd("~/Dropbox/Heise/U19-Ab/weight/chr1_microsatellite/")

pheno=read.csv("microsatellite_queries/preCTCTT_kmer_query1_1.csv")
score=read.csv("chr1_38-41_haplotypes.csv")
score=score[-2]

#set number of ID columns in phenotype sheet and QTL chromosome
hc=2
chromosome=1

num=length(pheno)-hc
# num=c(1:2,11:16,49:50)

pheno=merge(pheno,score)

colnames(pheno)[length(pheno)]=paste0("chr",chromosome,".status")

# #add mx1 status
# mx1=read.csv("~/Dropbox/Heise/CC/mx1_status.csv")
# mx1=mx1[c(1,10)]
# pheno=merge(pheno,mx1)
# colnames(pheno)[length(pheno)]="mx1.status"

# graph to view phenoa distribution #
abs=colnames(pheno)[(hc+1):length(pheno)]
par(mfrow = c(2,4))

for(i in 1:num)
{
  response.var<-pheno[,i+hc]
  hist(response.var, main=abs[i])
}

#phenotype distribution
for(i in 1:num)
{
  response.var=pheno[c(1:hc,i+hc)]

  g=ggplot(pheno, aes_string(x=abs[i],y=paste0("reorder(strain,",abs[i],")")))
  print(
  g+geom_point(aes_string(color=paste0("chr",chromosome,".status")))+theme_minimal()+
      labs(title=paste(abs[i]),y="strain",x="response variable")
  # +theme(legend.position='none')
  )
}

#phenotype vs. genotype
inheritance.model=paste0("chr",chromosome,".status")

for(i in num)
{
response.var=abs[i]

g=ggplot(pheno, aes_string(inheritance.model,response.var))
print(
g+labs(x=paste("Allele Score at QTL under",inheritance.model,"model",sep=" "),y=paste(response.var))+
  geom_jitter(width=0.1,aes(colour = strain))+theme_minimal()+theme(legend.position='none')
)
}
