library(doBy)
library(dplyr)
library(ggplot2)

setwd("~/Dropbox/Heise/U19-Ab/")

#load in data with phenotype
pheno=read.csv("weight/peak_weight_loss.csv")
mice="RIX" #or "line" or "strain" or whatever your column is called

#get allele scores (and mx1 status)
score=get.haplos2(chrom,start,end)
locus.name="chr1"

pheno=add.stat.gene(pheno,score,locus.name) #will duplicate individuals if their strain has a recombination/is heterozygous and you use gethaplos instead of gethaplos2
pheno=add.stat.mx1(pheno)

#get table of pairwise founder combinations (for RIX data)
founder.table=table(pheno[c(paste0("dam.founder_",locus.name),paste0("sire.founder_",locus.name))])
founder.table
# write.table(founder.table,"rixtable.csv",sep=",",col.names=NA)

#plot founder alleles vs. phenotype (really only useful for RI data)
inheritance.model="dam.founder_chr1"
response.var="lowest.weight"

founder.colors <- c(A="#F0F000", B="#808080", C="#F08080", D="#1010F0", E="#00A0F0", F="#00A000", G="#F00000", H="#9000E0")
founder.names <- c("AJ", "C57BL6J", "129S1", "NOD", "NZO", "CAST", "PWK", "WSB")

g=ggplot(pheno, aes_string(inheritance.model,response.var))
g+labs(x=paste(inheritance.model),y=paste(response.var))+
  geom_boxplot(aes(fill=get(inheritance.model)))+scale_fill_manual(values=founder.colors)+
  theme_minimal()+theme(legend.position='none')

#plot score based on haplotype vs. phenotype
inheritance.model="combo_chr1"
response.var="lowest.weight"
facet.var="combo_mx1"

g=ggplot(pheno, aes_string(inheritance.model,response.var))

#basic plot
g+labs(x=paste("allele score with",inheritance.model,"model",sep=" "),y=paste(response.var))+
  geom_jitter(width=0.2,aes_string(colour = inheritance.model))+theme_minimal()+theme(legend.position='none')

#faceted by mx1 status
g+labs(x=paste("Allele Score under",inheritance.model,"model",sep=" "),y=paste(response.var))+
  geom_jitter(width=0.2,aes(colour = get(mice)),na.rm=T)+
  # ylim(60,110)+
  theme_minimal() + guides(color=guide_legend(title="Mx1 allele combination")) +
  facet_wrap(facet.var) + theme(legend.position = "none") +
  ggtitle(paste(inheritance.model,"on",response.var,"grouped by",facet.var,sep=" "))

#to look at one haplotype subset only
haplotype="0_0"
locus.call="combo_mx1"

pheno.sub=subset(pheno,pheno[,names(pheno)==locus.call]==haplotype)
ggplot(pheno.sub, aes_string(inheritance.model, response.var))+
  labs(x=paste("Allele Score at QTL under",inheritance.model,"model",sep=" "),y=paste(response.var))+
  geom_jitter(width=0.3,aes(colour = get(mice)))+
  # xlim(c("0","1","2","3","4"))+
  theme_minimal()+theme(legend.position='none')+
  ggtitle(paste("combo_mx1",haplotype,"only",sep=" "))

######## strain averages
preserve.scores=c("combo_chr1","combo_mx1")
preserve.vars=c("lowest.weight")

preserve.sco.form = paste(preserve.scores,collapse="+")
preserve.var.form = paste(preserve.vars,collapse="+")
avg.formula=paste0(preserve.var.form,"~",mice,"+",preserve.sco.form)

pheno.avg=summaryBy(as.formula(avg.formula),data=pheno,FUN=mean,na.rm=T)
colnames(pheno.avg)[(length(pheno.avg)-length(preserve.vars)+1):length(pheno.avg)]=preserve.vars


g=ggplot(pheno.avg, aes_string(inheritance.model,response.var))


g+labs(x=paste("Allele Score under",inheritance.model,"model",sep=" "),y=paste(response.var))+
  geom_jitter(width=0.2,aes_string(colour = preserve.scores[1]))+theme_minimal()+
  ggtitle(paste("strain averages:",inheritance.model,"on",response.var,"grouped by",facet.var,sep=" ")) +
  guides(color=guide_legend(title="Mx1 allele combination")) +
  facet_wrap(facet.var) + theme(legend.position = "none")

## mx1 subset only
haplotype="0_0"
locus.call="combo_mx1"

pheno.sub=subset(pheno.avg,pheno.avg[,names(pheno.avg)==locus.call]==haplotype)

g+labs(x=paste("Allele Score under",inheritance.model,"model",sep=" "),y=paste(response.var))+
  geom_jitter(width=0.2,aes_string(colour = preserve.scores[1]))+theme_minimal()+
  # theme(legend.position='none')+
  ggtitle(paste("strain averages:",inheritance.model,"on",response.var,"(",locus.call,"=",haplotype,"strains only",")",sep=" "))

