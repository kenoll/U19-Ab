library(pwr)

#follow up from locus_allele_effect.R
inheritance.model="additive"
response.var="lowest.weight"

#all strains
g=ggplot(pheno.avg, aes_string(inheritance.model,response.var))
g+labs(x=paste("Allele Score at QTL under",inheritance.model,"model",sep=" "),y=paste(response.var))+
  geom_jitter(width=0.2,aes_string(colour = inheritance.model))+theme_minimal()+theme(legend.position='none')

#subset of mx1 dom/dom strains only
haplotype="0"
pheno.avg.mx1=subset(pheno.avg,pheno.avg$sum.mx1==haplotype)
g=ggplot(pheno.avg.mx1, aes_string(inheritance.model,response.var))
g+labs(x=paste("Allele Score at QTL under",inheritance.model,"model",sep=" "),y=paste(response.var))+
  geom_jitter(width=0.2,aes_string(colour = inheritance.model))+theme_minimal()+theme(legend.position='none')+
  ggtitle(paste("Weight loss based on Chr1 QTL allele Mx1",haplotype,"strains only",sep=" "))

#get means and sd for power calculations on G*power
score1=0
mean(subset(pheno.avg.mx1$lowest.weight,pheno.avg.mx1$additive==score1))
sd(subset(pheno.avg.mx1$lowest.weight,pheno.avg.mx1$additive==score1))

score2=3
mean(subset(pheno.avg.mx1$lowest.weight,pheno.avg.mx1$additive==score2))
sd(subset(pheno.avg.mx1$lowest.weight,pheno.avg.mx1$additive==score2))

# #power calcs in R
# pwr.t.test

#linear models comparison
plot(lowest.weight~additive,pheno.avg.mx1)
reg=lm(lowest.weight~additive,pheno.avg.mx1)
summary(reg)

pheno.mx1=subset(pheno,pheno$sum.mx1=="0")
library(lme4)
model_base<-lmer(lowest.weight~(1|RIX),data=pheno.mx1)
model_full<-lmer(lowest.weight~(1|RIX)+additive, data=pheno.mx1)
anova(model_base, model_full)

######
dat=pheno
library(dplyr)

# strain.mean=aggregate(lowest.weight~RIX,data=dat,FUN=mean)
#   colnames(strain.mean)[2]="mean"
strain.n=aggregate(lowest.weight~RIX,data=dat,FUN=length)
  colnames(strain.n)[2]="n"
strain.sd=aggregate(lowest.weight~RIX,data=dat,FUN=sd)
  colnames(strain.sd)[2]="sd"
# strain.var=merge(strain.mean,strain.sd)
# strain.var=merge(strain.var,strain.n)
# strain.var = strain.mean %>% merge(strain.sd) %>% merge(strain.n) %>% merge(pheno.avg)
se <- function(x) sd(x)/sqrt(length(x))
strain.se=aggregate(lowest.weight~RIX,data=dat,FUN=se)
  colnames(strain.se)[2]="se"
strain.var = strain.n %>% merge(strain.sd) %>% merge(strain.se) %>% merge(pheno.avg)

strain.var$additive = strain.var$additive %>% as.character %>% as.factor

library(forcats)

strain.var %>% 
  mutate(ordering = as.numeric(additive) + lowest.weight*.001,
         RIX = fct_reorder(RIX, ordering, .desc = T)) %>% 
  ggplot(aes(x=RIX, y=lowest.weight, color=additive)) + geom_point() +
    geom_errorbar(aes(ymin=lowest.weight-se, ymax=lowest.weight+se), width=.1) +
    geom_point() + facet_wrap("sum.mx1",scales="free_x") + theme_minimal() +
    theme(axis.text.x=element_text(angle=60,vjust=1,hjust=1))

strain.var %>% 
  mutate(ordering = as.numeric(additive) + lowest.weight*.001,
         RIX = fct_reorder(RIX, ordering, .desc = T)) %>% 
  subset(strain.var$sum.mx=="0") %>%
  ggplot(aes(x=RIX, y=lowest.weight, color=additive)) + geom_point() +
    geom_errorbar(aes(ymin=lowest.weight-se, ymax=lowest.weight+se), width=.1) +
    geom_point() + theme_minimal() +
    theme(axis.text.x=element_text(angle=60,vjust=1,hjust=1))


pheno.val = pheno.avg %>% subset(additive == "0" | additive == "3") %>% subset(sum.mx1=="0")
  pheno.val$RIX=as.character(pheno.val$RIX)
  RIX_sep <- data.frame(do.call("rbind", strsplit(pheno.val$RIX,"x")))
  colnames(RIX_sep)[1:2]=c("dam","sire")
  pheno.val=cbind(RIX_sep,pheno.val)
  
  write.csv(pheno.val,"qtl_info/chr1_validation_RIXs.csv",row.names=F)
