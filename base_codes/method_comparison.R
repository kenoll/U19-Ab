setwd("~/Dropbox/Heise/U19-Ab/antibody")
library(gridExtra)

methods = merge(aucdata,lastpos) %>% merge(halfmax)
methods = methods[c(1:3,6,13:15)]

abs=c("IgG1" ,  "IgG2ac" , "IgG2b" , "IgG3" ,  "IgM"  ,  "TotalG")
days=c("7","10","15","45")

for (k in 1:length(abs))
{
  iso=abs[k]
  dat=subset(methods,methods$isotype==iso)
 
  png(file.path(paste0("comparison/","method_histo_",iso)),width=700,height=700)
  par(mfrow = c(4,3))
      for(i in 1:length(days))
      {
        phenotype=subset(dat,dat$day==days[i])
        hist(phenotype$AUC, main=paste("day",days[i],iso,"AUC",sep=" "), col="papayawhip",xlab="")
        hist(phenotype$last_positive, main=paste("day",days[i],iso,"lastpos",sep=" "), col="lemonchiffon",xlab="")
        hist(phenotype$halfmax, main=paste("day",days[i],iso,"halfmax",sep=" "), col="mintcream",xlab="")
      }
  dev.off()
}




her=read.csv("comparison/heritability_estimates_methods.csv")
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

png(file.path(paste0("comparison/","method_heritability_plots")),width=1200,height=800)
grid.arrange(grobs = plot.list, ncol=3)
dev.off()