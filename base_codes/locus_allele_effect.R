setwd("~/Dropbox/Heise/U19-Ab/")

#load in data with phenotype
pheno=read.csv("weight/peak_weight_loss.csv")

#load in allele scores (see locus_haplotype.R)
locus="chr1_38-41_scores"
score=read.csv(paste0("qtl_info/",locus,".csv"))

# code with  score allele (e.g. for graphing)
pheno$RIX=as.character(pheno$RIX)

RIX_sep <- data.frame(do.call("rbind", strsplit(pheno$RIX,"x")))
colnames(RIX_sep)[1:2]=c("dam","sire")
pheno=cbind(RIX_sep,pheno)

#if CC names are in old alias, can use c(2,3)
score=score[c(1,3)]

colnames(score)[1]="sire"
pheno=merge(pheno,score)
colnames(pheno)[length(pheno)]="sire.score"

colnames(score)[1]="dam"
pheno=merge(pheno,score)
colnames(pheno)[length(pheno)]="dam.score"

pheno$additive=rowSums(pheno[c("dam.score","sire.score")])
pheno$additive=as.character(pheno$additive)
pheno$additive=as.factor(pheno$additive)

#plot
inheritance.model="additive"
response.var="lowest.weight"

g=ggplot(pheno, aes_string(inheritance.model,response.var))
g+labs(x=paste("Allele Score at QTL under",inheritance.model,"model",sep=" "),y=paste(response.var))+
  geom_jitter(width=0.2,aes_string(colour = inheritance.model))+theme_minimal()+theme(legend.position='none')

# ggsave(paste0("qtls/",locus,"_",response.var,"_",inheritance.model,"_plot.png"),width=5,height=3)

# g+labs(x=paste("Allele Score at QTL under",inheritance.model,"model",sep=" "),y=paste(response.var))+
#   geom_boxplot(aes_string(colour = inheritance.model))+theme_minimal()+theme(legend.position='none')

######################################################

      ###### to drop any ambiguous calls (heterozygosity or recombination within QTL)
      pheno=read.csv("peak_weight_loss.csv")
      pheno$RIX=as.character(pheno$RIX)
      RIX_sep <- data.frame(do.call("rbind", strsplit(pheno$RIX,"x")))
      colnames(RIX_sep)[1:2]=c("dam","sire")
      pheno=cbind(RIX_sep,pheno)
      
      #load in allele scores (see locus_haplotype.R)
      locus="chr1_38-41"
      score=read.csv(paste0("qtl_info/",locus,"_norecomb_haplotypes.csv"))
      
      #if CC names are in old alias, can use -1
      score=score[-2]
      
      colnames(score)[1]="sire"
      pheno=merge(pheno,score)
      colnames(pheno)[length(pheno)]="sire.score"
      colnames(pheno)[(length(pheno)-1)]="sire.founder"
      
      colnames(score)[1]="dam"
      pheno=merge(pheno,score)
      colnames(pheno)[length(pheno)]="dam.score"
      colnames(pheno)[(length(pheno)-1)]="dam.founder"
      
      pheno$additive=rowSums(pheno[c("dam.score","sire.score")])
      pheno$additive=as.character(pheno$additive)
      pheno$additive=as.factor(pheno$additive)
      
      pheno$combo.score=paste(pheno$dam.score,pheno$sire.score,sep="_")
      pheno$combo.score=gsub("0_1","1_0",pheno$combo.score)
      pheno$combo.score=gsub("0_2","2_0",pheno$combo.score)
      pheno$combo.score=gsub("1_2","2_1",pheno$combo.score)
      
      #plot
      inheritance.model="combo.score"
      response.var="lowest.weight"
      
      g.sub=ggplot(pheno, aes_string(inheritance.model,response.var))
      
      g.sub+labs(x=paste("Allele Score at QTL under",inheritance.model,"model",sep=" "),y=paste(response.var))+
        geom_jitter(width=0.2,aes_string(colour = inheritance.model))+
        theme_minimal()+theme(legend.position='none')
      
      # ggsave(paste0("qtls/",locus,"_",response.var,"_",inheritance.model,"_plot_norecombinations.png"),width=5,height=3)
      
      # g.sub+labs(x=paste("Allele Score at QTL under",inheritance.model,"model",sep=" "),y=paste(response.var))+
      #   geom_jitter(width=0.2,aes(colour = RIX))+
      #   # geom_jitter(width=0.2,aes(colour = sire.founder))+
      #   theme_minimal()+theme(legend.position='none')


######## add mx1 #########
mx1=read.csv("~/Dropbox/Heise/CC/mx1_status.csv")
mx1=mx1[c(1,10)]

colnames(mx1)[1]="dam"
pheno=merge(pheno,mx1)
colnames(pheno)[length(pheno)]="dam.mx1"

colnames(mx1)[1]="sire"
pheno=merge(pheno,mx1)
colnames(pheno)[length(pheno)]="sire.mx1"

pheno$dam.mx1=as.character(pheno$dam.mx1)
pheno$sire.mx1=as.character(pheno$sire.mx1)
pheno$combo.mx1=paste(pheno$dam.mx1,pheno$sire.mx1,sep="_")
pheno$combo.mx1=gsub("0_1","1_0",pheno$combo.mx1)
pheno$combo.mx1=gsub("0_0.5","0.5_0",pheno$combo.mx1)

pheno$dam.mx1=as.numeric(pheno$dam.mx1)
pheno$sire.mx1=as.numeric(pheno$sire.mx1)
pheno$sum.mx1=rowSums(pheno[c("dam.mx1","sire.mx1")])

pheno$sum.mx1=as.factor(pheno$sum.mx1)
pheno$dam.mx1=as.factor(pheno$dam.mx1)
pheno$sire.mx1=as.factor(pheno$sire.mx1)


inheritance.model="additive"
response.var="lowest.weight"

g.sub.mx1=ggplot(pheno, aes_string(inheritance.model,response.var))

#mx1 allele score combo: parent of origin ordered (0/1 different than 1/0)
g.sub.mx1+labs(x=paste("Allele Score at QTL under",inheritance.model,"model",sep=" "),y=paste(response.var))+
  geom_jitter(width=0.3,aes(colour = interaction(dam.mx1,sire.mx1,sep=" / ",lex.order = T)))+
  theme_minimal() + guides(color=guide_legend(title="Mx1 allele combination"))

#mx1 allele score combo: dam vs sire unordered (0/1 same as 1/0)
g.sub.mx1+labs(x=paste("Allele Score at QTL under",inheritance.model,"model",sep=" "),y=paste(response.var))+
  geom_jitter(width=0.3,size=1,aes(color = combo.mx1),na.rm=T)+
  theme_minimal() + guides(color=guide_legend(title="Mx1 allele combination"))

#faceted by mx1 status
g.sub.mx1+labs(x=paste("Allele Score at QTL under",inheritance.model,"model",sep=" "),y=paste(response.var))+
  geom_jitter(width=0.3,aes(colour = RIX),na.rm=T)+ylim(60,110)+
  theme_minimal() + guides(color=guide_legend(title="Mx1 allele combination")) +
  facet_wrap("combo.mx1") + theme(legend.position = "none") +
  ggtitle("Chr1 QTL effects, grouped by Mx1 allele combination")

#to look at one mx1 subset only
haplotype="0_0"

pheno.mx1=subset(pheno,pheno$combo.mx1==haplotype)
ggplot(pheno.mx1, aes(additive, lowest.weight))+
  labs(x=paste("Allele Score at QTL under",inheritance.model,"model",sep=" "),y=paste(response.var))+
  geom_jitter(width=0.3,aes(colour = RIX))+xlim(c("0","1","2","3","4"))+
  theme_minimal()+theme(legend.position='none')+
  ggtitle(paste("Mx1",haplotype,"only",sep=" "))

######## strain averages
library(doBy)
pheno.avg=summaryBy(lowest.weight~RIX+additive+sum.mx1, data=pheno, FUN=mean, na.rm=T)
colnames(pheno.avg)[length(pheno.avg)]="lowest.weight"

g=ggplot(pheno.avg, aes_string(inheritance.model,response.var))
g+labs(x=paste("Allele Score at QTL under",inheritance.model,"model",sep=" "),y=paste(response.var))+
  geom_jitter(width=0.2,aes_string(colour = "sum.mx1"))+theme_minimal()+geom_text(aes(label=RIX))

  # theme(legend.position='none')+
  ggtitle("Weight loss based on Chr1 QTL allele Mx1, strain averages")

      
      ## mx1 subset only
      haplotype="0"
      
      pheno.avg.mx1=subset(pheno.avg,pheno.avg$sum.mx1==haplotype)
      g=ggplot(pheno.avg.mx1, aes_string(inheritance.model,response.var))
      g+labs(x=paste("Allele Score at QTL under",inheritance.model,"model",sep=" "),y=paste(response.var))+
        geom_jitter(width=0.2,aes_string(colour = inheritance.model))+theme_minimal()+theme(legend.position='none')+
        ggtitle(paste("Weight loss based on Chr1 QTL allele Mx1",haplotype,"strains only",sep=" "))

