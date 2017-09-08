setwd("~/Dropbox/Heise/U19-Ab/antibody/")
pheno=read.csv("data_bin/d10/dat.10.auc.csv")

pheno = pheno %>% add.stat.mx1
preserve.scores=c("combo_mx1")
preserve.vars=c("TotalG","IgG1","IgG2ac","IgG2b","IgG3","IgM")

preserve.sco.form = paste(preserve.scores,collapse="+")
preserve.var.form = paste(preserve.vars,collapse="+")
avg.formula=paste0(preserve.var.form,"~",mice,"+",preserve.sco.form)
pheno.avg=summaryBy(as.formula(avg.formula),data=pheno,FUN=mean,na.rm=T)
colnames(pheno.avg)[(length(pheno.avg)-length(preserve.vars)+1):length(pheno.avg)]=preserve.vars

inheritance.model="combo_mx1"
response.var="TotalG"

g=ggplot(pheno.avg, aes_string(inheritance.model,response.var))

g+labs(x=paste("Allele Score under",inheritance.model,"model",sep=" "),y=paste(response.var))+
  geom_jitter(width=0.2,aes_string(colour = preserve.scores[1]),na.rm=T)+theme_minimal()+
  # theme(legend.position='none')+
  ggtitle(paste("strain averages:",inheritance.model,"on",response.var,sep=" "))

strains=subset(pheno.avg,pheno.avg$TotalG>5 | pheno.avg$TotalG<2.5)
strains=subset(strains,strains$combo_mx1 %in% c("0_0","0_1"))

g=ggplot(strains, aes_string(inheritance.model,response.var))
g+labs(x=paste("Allele Score under",inheritance.model,"model",sep=" "),y=paste(response.var))+
  geom_jitter(width=0.2,aes_string(colour = preserve.scores[1]),na.rm=T)+theme_minimal()+
  ggtitle(paste("strain averages:",inheritance.model,"on",response.var,sep=" "))


#
rixs=read.csv("~/Dropbox/Heise/CC/U19_RIX_list.csv")
rixs = rixs %>% add.stat.mx1
table(rixs$combo_mx1)

strains.old=merge(rixs[3:4],strains[1:3])
