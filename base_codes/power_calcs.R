library(pwr)
library(doBy)
library(dplyr)
library(ggplot2)
library(lme4)

#follow up from locus_allele_effect.R
setwd("~/Dropbox/Heise/U19-Ab/")

pheno=read.csv("weight/weights_pers_2016_10.csv")
pheno=pheno[c(1:6,13)]
pheno$Day=as.character(pheno$Day)
pheno$Day=as.numeric(pheno$Day)
pheno=subset(pheno,pheno$Day >= 7)
pheno=pheno[complete.cases(pheno),]

# pheno=read.csv("weight/peak_weight_loss.csv")

score=read.csv("qtl_info/chr1_38-41_scores.csv")
pheno = pheno %>% add.stat.gene(score,"chr1") %>% add.stat.mx1

pheno.avg=summaryBy(D7per~RIX+add_chr1+add_mx1, data=pheno, FUN=mean, na.rm=T)
colnames(pheno.avg)[length(pheno.avg)]="D7per"

inheritance.model="add_chr1"
response.var="D7per"


#all strains
g=ggplot(pheno.avg, aes_string(inheritance.model,response.var))
g+labs(x=paste("Allele Score at QTL under",inheritance.model,"model",sep=" "),y=paste(response.var))+
  geom_jitter(width=0.2,aes_string(colour = inheritance.model))+theme_minimal()+theme(legend.position='none')

#subset of mx1 dom/dom strains only
haplotype="0"
pheno.avg.mx1=subset(pheno.avg,pheno.avg$add_mx1==haplotype)
g=ggplot(pheno.avg.mx1, aes_string(inheritance.model,response.var))
g+labs(x=paste("Allele Score at QTL under",inheritance.model,"model",sep=" "),y=paste(response.var))+
  geom_jitter(width=0.2,aes_string(colour = inheritance.model))+theme_minimal()+theme(legend.position='none')+
  ggtitle(paste("Weight loss based on Chr1 QTL allele Mx1",haplotype,"strains only",sep=" "))

#get means and sd for power calculations on G*power
score1=0
mean(subset(pheno.avg.mx1$D7per,pheno.avg.mx1$add_chr1==score1))
sd(subset(pheno.avg.mx1$D7per,pheno.avg.mx1$add_chr1==score1))

score2=3
mean(subset(pheno.avg.mx1$D7per,pheno.avg.mx1$add_chr1==score2))
sd(subset(pheno.avg.mx1$D7per,pheno.avg.mx1$add_chr1==score2))

# #power calcs in R
# pwr.t.test

#linear models comparison
plot(D7per~add_chr1,pheno.avg.mx1)
reg=lm(D7per~add_chr1,pheno.avg.mx1)
summary(reg)

pheno.mx1=subset(pheno,pheno$add_mx1=="0")
model_base<-lmer(D7per~(1|RIX),data=pheno.mx1)
model_full<-lmer(D7per~(1|RIX)+add_chr1, data=pheno.mx1)
anova(model_base, model_full)

######
dat=pheno

# strain.mean=aggregate(D7per~RIX,data=dat,FUN=mean)
#   colnames(strain.mean)[2]="mean"
strain.n=aggregate(D7per~RIX,data=dat,FUN=length)
  colnames(strain.n)[2]="n"
strain.sd=aggregate(D7per~RIX,data=dat,FUN=sd)
  colnames(strain.sd)[2]="sd"
# strain.var=merge(strain.mean,strain.sd)
# strain.var=merge(strain.var,strain.n)
# strain.var = strain.mean %>% merge(strain.sd) %>% merge(strain.n) %>% merge(pheno.avg)
se <- function(x) sd(x)/sqrt(length(x))
strain.se=aggregate(D7per~RIX,data=dat,FUN=se)
  colnames(strain.se)[2]="se"
strain.var = strain.n %>% merge(strain.sd) %>% merge(strain.se) %>% merge(pheno.avg)

strain.var$add_chr1 = strain.var$add_chr1 %>% as.character %>% as.factor

strain.var.sub = strain.var %>% subset(add_chr1 == "0" | add_chr1 == "3") %>% subset(add_mx1=="0")
strain.var.3 = strain.var %>% subset(add_chr1 == "3") %>% subset(add_mx1=="0")
strain.var.0 = strain.var %>% subset(add_chr1 == "0") %>% subset(add_mx1=="0")
strain.var.3
strain.var.0

dodge=position_dodge(width=1)
eb="sd"

ggplot(strain.var.sub, aes_string(inheritance.model,response.var,group=response.var,color=inheritance.model))+
  labs(x=paste("Allele Score at QTL under",inheritance.model,"model",sep=" "),y=paste(response.var))+
  theme_minimal()+theme(legend.position='none')+
  geom_point(position=dodge,aes(size=n))+
  geom_errorbar(aes(ymin=get(response.var)-get(eb), ymax=get(response.var)+get(eb)), position=dodge,width=.1)


###############

# library(forcats)
# 
# strain.var %>% 
#   mutate(ordering = as.numeric(add_chr1) + D7per*.001,
#          RIX = fct_reorder(RIX, ordering, .desc = T)) %>% 
#   ggplot(aes(x=RIX, y=D7per, color=add_chr1)) + geom_point() +
#     geom_errorbar(aes(ymin=D7per-se, ymax=D7per+se), width=.1) +
#     geom_point() + facet_wrap("add_mx1",scales="free_x") + theme_minimal() +
#     theme(axis.text.x=element_text(angle=60,vjust=1,hjust=1))
# 
# strain.var %>% 
#   mutate(ordering = as.numeric(add_chr1) + D7per*.001,
#          RIX = fct_reorder(RIX, ordering, .desc = T)) %>% 
#   subset(strain.var$sum.mx=="0") %>%
#   ggplot(aes(x=RIX, y=D7per, color=add_chr1)) + geom_point() +
#     geom_errorbar(aes(ymin=D7per-se, ymax=D7per+se), width=.1) +
#     geom_point() + theme_minimal() +
#     theme(axis.text.x=element_text(angle=60,vjust=1,hjust=1))
# 
# 
# pheno.val = pheno.avg %>% subset(add_chr1 == "0" | add_chr1 == "3") %>% subset(add_mx1=="0")
#   pheno.val$RIX=as.character(pheno.val$RIX)
#   RIX_sep <- data.frame(do.call("rbind", strsplit(pheno.val$RIX,"x")))
#   colnames(RIX_sep)[1:2]=c("dam","sire")
#   pheno.val=cbind(RIX_sep,pheno.val)
#   
#   write.csv(pheno.val,"qtl_info/chr1_validation_RIXs.csv",row.names=F)
#   

