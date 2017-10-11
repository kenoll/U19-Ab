setwd("~/Dropbox/Heise/U19-Ab/antibody")
library(gridExtra)
library(ggplot2)
library(reshape2)
library(dplyr)

#### PHENOTYPE DISTRIBUTION ####

aucdata=read.csv("data_bin/auc_transformed_alldays.csv")
  aucdata=aucdata[c(1,4,5,11:16)]
  aucdata=melt(aucdata,id.vars=c("RIX","ID","day"),variable.name="isotype",value.name="AUC")
lastpos=read.csv("data_bin/lastpos_transformed_alldays.csv")
  lastpos=lastpos[c(1,4,5,11:16)]
  lastpos=melt(lastpos,id.vars=c("RIX","ID","day"),variable.name="isotype",value.name="lastpos")
halfmax=read.csv("data_bin/halfmax_transformed_alldays.csv")
  halfmax=halfmax[c(1,4,5,11:16)]
  halfmax=melt(halfmax,id.vars=c("RIX","ID","day"),variable.name="isotype",value.name="halfmax")

methods = merge(aucdata,lastpos) %>% merge(halfmax)

abs=c("IgG1" ,  "IgG2ac" , "IgG2b" , "IgG3" ,  "IgM"  ,  "TotalG")
days=c("7","10","15","45")

for (k in 1:length(abs))
{
  iso=abs[k]
  dat=subset(methods,methods$isotype==iso)
 
  png(file.path(paste0("comparison/method_histo_",iso,".png")),width=600,height=600)
  par(mfrow = c(4,3))
      for(i in 1:length(days))
      {
        phenotype=subset(dat,dat$day==days[i])
        hist(phenotype$AUC, main=paste("day",days[i],iso,"AUC",sep=" "), col="papayawhip",xlab="")
        hist(phenotype$lastpos, main=paste("day",days[i],iso,"lastpos",sep=" "), col="lemonchiffon",xlab="")
        hist(phenotype$halfmax, main=paste("day",days[i],iso,"halfmax",sep=" "), col="mintcream",xlab="")
      }
  dev.off()
}












#### HERITABILITY ####

her.auc=read.csv("comparison/heritability_estimates_auc_alldays.csv")
her.halfmax=read.csv("comparison/heritability_estimates_halfmax_alldays.csv")
her.lastpos=read.csv("comparison/heritability_estimates_lastpos_alldays.csv")
her=rbind(her.auc,her.halfmax,her.lastpos)

her$day=as.factor(her$day)

plot.list=NULL
for (k in 1:length(abs))
{
  iso=abs[k]
  dat=subset(her,her$isotype==iso)

  plot.list[[k]]=
  ggplot(dat,aes(x=day,y=(upper+lower)/2,fill=method),group=method)+
  geom_col(position=position_dodge())+
  # geom_errorbar(aes(ymin=lower, ymax=upper),position=position_dodge())+
  scale_fill_brewer(palette="Set2")+ggtitle(paste0(iso," heritability estimates"))
}

png("comparison/method_heritability_plots.png",width=1200,height=800)
grid.arrange(grobs = plot.list, ncol=3)
dev.off()
