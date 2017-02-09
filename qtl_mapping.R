#### setup: load libraries, data format, etc ####
setwd("~/Dropbox/Heise/ELISA Antibody/")

#set date
run.date="2017_01_30"
date.to.load="2016_01_16"

#load libraries 
library(ggplot2)

#load in transformed data (see "Antibody data for mapping pipeline.R")
dat.7=read.csv("~/Dropbox/Heise/ELISA Antibody/qtls/d7/dat.7.csv")
dat.10=read.csv("~/Dropbox/Heise/ELISA Antibody/qtls/d10/dat.10.csv")
dat.15=read.csv("~/Dropbox/Heise/ELISA Antibody/qtls/d15/dat.15.csv")
dat.45=read.csv("~/Dropbox/Heise/ELISA Antibody/qtls/d45/dat.45.csv")

#set list of anitbody isotypes in the column order they appear in
abs=c("IgG1" ,  "IgG2ac" , "IgG2b" , "IgG3" ,  "IgM"  ,  "TotalG")

col.list=c("RIX" , "ID" , "day" , "assay_date", "IgG1" ,"IgG2ac","IgG2b" ,    
           "IgG3","IgM","TotalG")

#### lists ####
day.list=list(7,10,15,45)
dat.list=list(dat.7,dat.10,dat.15,dat.45)

#### phenotypic distribution plots ####

for(k in 1:4)
{
  day=day.list[[k]]
  dat=dat.list[[k]]
  
  for(i in 1:6)
  {
    phenotype=dat[c(1:4,i+4)]
    phenotype=phenotype[complete.cases(phenotype),]  
    
    g<-ggplot(phenotype, aes_string(x=abs[i],y=paste0("reorder(RIX,",abs[i],")")))
    g+geom_point(aes(colour=assay_date))+theme_classic()+
      labs(title=paste("Day",day.list[[k]],abs[i],sep=" "),y="RIX",x="AUC")
    
    ggsave(file.path(paste("qtls/d",day,sep=""),paste("b6_pheno_dist_batch_",abs[i],"_d",day.list[[k]],".jpg",sep="")), width=6, height=9)
  }
}


#geom_point(aes(colour=assay_date))
#+ geom_point(colour="purple")

# # test individually # #
# dat=dat.10
# i=1
# 
# phenotype=dat[c(1:4,i+4)]
# phenotype=phenotype[complete.cases(phenotype),]  
# 
# g<-ggplot(phenotype, aes_string(x=abs[i],y=paste0("reorder(RIX,",abs[i],")")))
# g+geom_point(aes(color=assay_date))+theme_classic()+
#   labs(title=paste("Day",day.list[[k]],abs[i],sep=" "),y="RIX",x="AUC")
#   # +  geom_text(aes(label=ID))

### phenotypic distribution plots / old code = loops for all isotypes of one day at a time ###
# 
# day=7
# # (paste0("dat.",day))
# dat=dat.7
# 
# for(i in 1:6)
# {
#   phenotype=dat[c(1:4,i+4)]
#   dat=dat[complete.cases(dat),]  
# 
# g<-ggplot(dat, aes_string(x=abs[i],y=paste0("reorder(RIX,",abs[i],")")))
# g+geom_point(aes(color=assay_date))+theme_classic()+
#   labs(title=paste("D10",abs[i],sep=" "),y="RIX",x="AUC")
# 
# ggsave(file.path(paste("qtls/d",day,sep=""),paste("pheno_dist_",abs[i],".jpg",sep="")), width=6, height=9)
# }
# 
# # "RIX"    "ID"     "day"    "virus"  "IgG1"   "IgG2ac" "IgG2b"  "IgG3"   "IgM"    "TotalG"


#### heritability calculations ####
heritability.bounds=NULL
N=3

for(k in 1:4)
{
  day=day.list[[k]]
  dat=dat.list[[k]]
  
  #set n as average number of replicates per strain
  
  for(i in 1:6)
  {
    a<-anova(lm(dat[,i+4]~dat$RIX))
    interclass.corr<-(a$"Mean Sq"[1] - a$"Mean Sq"[2])/(a$"Mean Sq"[1]+(N-1)*a$"Mean Sq"[2])
    Coef.Genet.Det<-(a$"Mean Sq"[1] - a$"Mean Sq"[2])/(a$"Mean Sq"[1]+((2*N)-1)*a$"Mean Sq"[2])
    h.b=c(interclass.corr,Coef.Genet.Det)
    heritability.bounds <- rbind(heritability.bounds, h.b)
  }
  
  rownames(heritability.bounds)=colnames(dat[5:10])
  heritability.bounds=as.data.frame(heritability.bounds)
  heritability.bounds[,3]=rep.int(day.list[[k]],nrow(heritability.bounds))
  colnames(heritability.bounds)=c("upper","lower","day")
  
  write.csv(heritability.bounds,
            file.path(paste("qtls/d",day,sep=""),paste("heritability_estimates_d",day.list[[k]],"_",run.date,".csv",sep="")),
  )
  heritability.bounds=NULL
}

#### correlation matrix ####
for(k in 1:4)
{
  day=day.list[[k]]
  dat=dat.list[[k]]
  
  V<-var(dat[,5:10], na.rm=T)
  C<-cov2cor(V)
  
  jpeg(file.path(paste("qtls/d",day,sep=""),paste("correlation_matrix_d",day.list[[k]],"_",run.date,".jpg",sep="")))
  heatmap(C, symm=TRUE, col = cm.colors(256))
  dev.off()
}

#### two-phenotype distributions ####
library(doBy)
dat=dat.10

Map.dat=summaryBy(IgG1+IgG2ac+IgG2b+IgG3+IgM+TotalG ~ 
                    RIX + day, data=dat, FUN=mean, na.rm=T)
colnames(Map.dat)[3:8] = gsub(".mean","",colnames(Map.dat[3:8]))
Map.dat$RIX=as.character(Map.dat$RIX)

#code with Mx1 status

RIX_sep <- data.frame(do.call("rbind", strsplit(Map.dat$RIX,"x")))
colnames(RIX_sep)[1:2]=c("dam","sire")

Map.dat=cbind(RIX_sep,Map.dat)

mx1=read.csv("~/Dropbox/Heise/CC/mx1_status.csv")
mx1=mx1[c(1,8)]

colnames(mx1)[1]="dam"
Map.dat=merge(Map.dat,mx1)
colnames(Map.dat)[11]="dam.mx1"

colnames(mx1)[1]="sire"
Map.dat=merge(Map.dat,mx1)
colnames(Map.dat)[12]="sire.mx1"

Map.dat=Map.dat[c(1,12,2,11,3:10)]

#graphs
Map.dat$RIX=as.factor(Map.dat$RIX)

phenotype=Map.dat[c(1:6,12,11)]
phenotype=phenotype[complete.cases(phenotype),]  

cbPalette <- c("royalblue","orchid","seagreen2")

g=ggplot(phenotype,aes(x=TotalG,y=IgM))
g+geom_point(size=4,aes(color=dam.mx1))+
  scale_color_manual(values=cbPalette)+
  geom_point(aes(color=sire.mx1))+
  theme_classic()+geom_text(aes(label=RIX),size=3,vjust=-1)+
  theme(legend.title=element_blank())


# rix.list=c("CC037xCC046","CC063xCC002","CC051xCC009","CC034xCC016",
#            "CC016xCC038","CC016xCC061","CC060xCC006")

# rix.list=("CC074xCC062")


rix.list=c("CC003xCC062")

rix.subset=Map.dat[which(Map.dat$RIX %in% rix.list),]

g=ggplot(phenotype,aes(x=TotalG,y=IgM))
g+geom_point(size=4,aes(color=dam.mx1))+
  scale_color_manual(values=cbPalette)+
  geom_point(aes(color=sire.mx1))+
  theme_classic()+geom_text(aes(label=RIX),size=3,vjust=-1)+
  theme(legend.title=element_blank())+
  geom_point(data=rix.subset,aes(x=TotalG,y=IgM),colour="red", size=5)


# ggsave("qtls/d45/totalg_v_igm_d7_highlighted.jpg",width=10,height=6)


#### mapping ####
#load in required libraries and data
library(doBy)
library(DOQTL)
load("~/Dropbox/Heise/ELISA Antibody/R codes and data/CCRIXb38F.Rdata")
# load(url("ftp://ftp.jax.org/MUGA/MM_snps.Rdata"))
load("~/Dropbox/Heise/ELISA Antibody/R codes and data/MM_snps.Rdata")

model.probs<-model.probs+(1e-20)
K=kinship.probs(model.probs)

## to graph QTL scans (no significance threshold yet)

for(k in c(1,4))
  # for(k in 1:4)
{
  day=day.list[[k]]
  dat=dat.list[[k]]
  
  Map.dat = summaryBy(IgG1+IgG2ac+IgG2b+IgG3+IgM+TotalG ~ 
                        RIX + day, data=dat, FUN=mean, na.rm=T)
  colnames(Map.dat)[3:8] = gsub(".mean","",colnames(Map.dat[3:8]))
  
  Map.dat$Sex="F"
  Map.dat$RIX=as.factor(Map.dat$RIX)
  
  row.names(Map.dat)<-Map.dat$RIX
  covar = data.frame(sex=as.numeric(Map.dat$Sex == "F"))
  rownames(covar)=rownames(Map.dat)
  
  png(file.path(paste("qtls/d",day,sep=""),paste("qtl_scan_d",day.list[[k]],"_",run.date,".png",sep="")),width=1000,height=500)
  par(mfrow = c(2,3))
  for(i in 1:6)
  {
    pheno=abs[i]
    qtl=scanone(pheno=Map.dat, pheno.col=pheno, addcovar=covar, probs=model.probs, K=K, snps=MM_snps)
    plot(qtl,main=paste("Day",day.list[[k]],pheno,sep=" "))
    saveRDS(qtl,file.path(paste("qtls/d",day,sep=""),paste("qtl_scan_d",day,"_",abs[i],"_",run.date,".Rdata",sep="")))
    #     save(qtl,file=file.path(paste("qtls/d",day,sep=""),paste("qtl_scan_d",day,"_",abs[i],".RData",sep="")))
  }
  dev.off()
}

# save PDFs w/ width of 11 and height of 6
# save PNGs w/ width of 1000 and height of 500

#### to graph with thresholds (takes a long time!) ####
#using a threshold of 0.9 for suggestive QTL

for(k in c(1,4))
  # for(k in 1:4)
{
  day=day.list[[k]]
  dat=dat.list[[k]]
  
  Map.dat = summaryBy(IgG1+IgG2ac+IgG2b+IgG3+IgM+TotalG ~ 
                        RIX + day, data=dat, FUN=mean, na.rm=T)
  colnames(Map.dat)[3:8] = gsub(".mean","",colnames(Map.dat[3:8]))
  
  Map.dat$Sex="F"
  Map.dat$RIX=as.factor(Map.dat$RIX)
  
  row.names(Map.dat)<-Map.dat$RIX
  covar = data.frame(sex=as.numeric(Map.dat$Sex == "F"))
  rownames(covar)=rownames(Map.dat)
  
  png(file.path(paste("qtls/d",day,sep=""),paste("qtl_scan_95_100_d",day.list[[k]],"_",run.date,".png",sep="")),width=1000,height=500)
  par(mfrow = c(2,3))
  for(i in 1:6)
  {
    pheno=abs[i]
    qtl=scanone(pheno=Map.dat, pheno.col=pheno, addcovar=covar, probs=model.probs, K=K, snps=MM_snps)
    
    perms = scanone.perm(pheno = Map.dat, pheno.col = abs[i], probs = model.probs, addcovar = covar, snps= MM_snps, nperm = 100)
    save(perms,file=file.path(paste("qtls/d",day,sep=""),paste("perms500_d",day.list[[k]],"_",abs[i],"_",run.date,".Rdata",sep="")))
    
    thr = quantile(perms, probs = 0.90)    
    
    plot(qtl,sig.thr=thr,main=paste("Day",day.list[[k]],pheno,sep=" "))
    
    
  }
  dev.off()
}












#### zoom in on chromosomes with significant/suggestive QTL to look at range and allele effects ####5
pheno="IgG"
qtl=scanone(pheno=Map.dat, pheno.col=pheno, addcovar=covar, probs=model.probs, K=K, snps=MM_snps)

chr=5

plot(qtl, main=paste("Day",day,pheno,sep=" "))
coefplot(qtl, chr=chr,legend=F)
coefplot_v2(qtl, chr=chr,main=paste("Day",day,pheno,"Chromosome",chr,sep=" "),legend=F)

#### mapping with saved qtl objects and permutations ####
# for(k in 1:4)
#{

date.to.load = "2017_01_16"

k=2

day=day.list[[k]]
dat=dat.list[[k]]

Map.dat = summaryBy(IgG1+IgG2ac+IgG2b+IgG3+IgM+TotalG ~ 
                      RIX + day, data=dat, FUN=mean, na.rm=T)
colnames(Map.dat)[3:8] = gsub(".mean","",colnames(Map.dat[3:8]))

Map.dat$Sex="F"
Map.dat$RIX=as.factor(Map.dat$RIX)

row.names(Map.dat)<-Map.dat$RIX
covar = data.frame(sex=as.numeric(Map.dat$Sex == "F"))
rownames(covar)=rownames(Map.dat)

png(file.path(paste("qtls/d",day,sep=""),paste("qtl_scan_95_d",day.list[[k]],"_",run.date,".png",sep="")),width=1000,height=500)
par(mfrow = c(2,3))
for(i in 1:6)
{
  pheno=abs[i]
  load(file=file.path(paste("qtls/d",day,sep=""),paste("qtl_scan_d",day.list[[k]],"_",abs[i],".Rdata",sep="")))
  
  #perms = 
  load(file=file.path(paste("qtls/d",day,sep=""),paste("perms500_d",day.list[[k]],"_",abs[i],"_",date.to.load,".Rdata",sep="")))
  
  thr = quantile(perms, probs = 0.95)    
  
  plot(qtl,sig.thr=thr,main=paste("Day",day.list[[k]],pheno,sep=" "))
}
dev.off()
}


### zoom in on chromosomes with significant/suggestive QTL to look at range and allele effects ###
day=7
pheno="IgG2ac"
load(file=file.path(paste("qtls/d",day,sep=""),paste("qtl_scan_d",day,"_",pheno,"_",date.to.load,".Rdata",sep="")))

chr=3

# plot(qtl, main=paste("Day",day,pheno,sep=" "))
# coefplot(qtl, chr=chr,legend=F)
coefplot_v2(qtl, chr=chr,main=paste("Day",day,pheno,"Chromosome",chr,sep=" "),legend=F)

# plot(lod[,3], coef[,colors[1,1]], type = "l", col = colors[1,3], lwd = 2,
#      ylim = c(-3, 3), xlab = 
#        paste("Chr", chr, "(Mb)"), ylab = "Founder Effects", axes = FALSE, ...)

# str(qtl$lod$A)

### to get regions of genome that fall under significance threshold
load(file=file.path(paste("qtls/d",day,sep=""),paste("perms100_d",day.list[[k]],"_",pheno,"_",run.date,".Rdata",sep="")))
thr = quantile(perms, probs = 0.95) 

#max(qtl$lod$A$lod)
qtl$lod$A[which.max(qtl$lod$A$lod),]
qtl$lod$A[(which.max(qtl$lod$A$lod)-50):(which.max(qtl$lod$A$lod)+20),]

#get confidence interval with 1.5 LOD drop
lod.peak=qtl$lod$A[which.max(qtl$lod$A$lod),"lod"]
lod.min=lod.peak-1.5
sig.range=subset(qtl$lod$A,qtl$lod$A$lod>lod.min & qtl$lod$A$chr==11)
sig.range[1,]
sig.range[nrow(sig.range),]

#compare to bayesian credible interval
bayesint(qtl,chr=11,prob=0.95,expandtomarkers=T)

#get haplotype probs for RIXs at the peaks - manually 
igg3.11.probs<-model.probs[,,qtl$lod$A[which.max(qtl$lod$A$lod),"marker"]]
write.csv(igg3.11.probs,"qtls/d10/igg3.11.qtl.probs.csv")


#to look at a peak other than the very max, input chromosome of interest
qtl.chr=subset(qtl$lod$A,qtl$lod$A$chr==5)
qtl.chr[which.max(qtl.chr$lod),]
qtl.chr[(which.max(qtl.chr$lod)-20):(which.max(qtl.chr$lod)+20),]

#to zoom in on qtl plot in region of interest
qtl.chr.sub=subset(qtl.chr,qtl.chr$pos>63 & qtl.chr$pos<75)
# plot(x=qtl.chr.sub$pos,y=qtl.chr.sub$lod)
ggplot(qtl.chr.sub, aes(pos, lod))+geom_point()+theme_minimal()


#11:70800000-72550000
#14 0.95 thr=6.707  / 0.90 thr= 6.481
#14 above .9 - 14:104754800-104996500

#### lms ####
#keyed mx1 status in with mx1.status.buildin code
#^^ requires Map.dat set up
str(dat.10)

lm.1=aov(TotalG~RIX,data=dat.10)
layout(matrix(c(1,2,3,4),2,2)) + plot(lm.1)
summary(lm.1)

lm.2=aov(TotalG~RIX+sum.mx1,data=dat.10)
layout(matrix(c(1,2,3,4),2,2)) + plot(lm.2)
summary(lm.2)

lm.3=lm(TotalG~sire.mx1+dam.mx1,data=dat.10)
# layout(matrix(c(1,2,3,4),2,2)) + plot(lm.3)
summary(lm.3)

lm.4=lm(TotalG~sum.mx1+sire.mx1+dam.mx1,data=dat.10)
# layout(matrix(c(1,2,3,4),2,2)) + plot(lm.3)
summary(lm.4)

lm.5=lm(TotalG~sum.mx1,data=dat.10)
layout(matrix(c(1,2,3,4),2,2)) + plot(lm.3)
summary(lm.5)


anova(lm.4,lm.3)

library(nlme)
m1.nlme = lme(TotalG~sum.mx1,random = ~ 1|RIX, data = dat.10,na.omit)

library(lme4)
m1.lme4 = lmer(TotalG~RIX + (1|sum.mx1),
               data = dat.10)


dat.10$sum.mx1=as.character(dat.10$sum.mx1)
dat.10$sum.mx1=as.factor(dat.10$sum.mx1)

g=ggplot(dat.10, aes(sum.mx1, IgG3))+labs(x="Allele Score",y="AUC")+
  geom_jitter(width=0.5,aes(colour = sum.mx1))+theme_minimal()+theme(legend.position='none')
g
g+geom_boxplot(aes(middle=mean(data), color=sum.mx1)) 
# g+geom_line(stat = "summary", fun.y = "mean", group="sum.mx1")


g=ggplot(dat.10, aes(sire.mx1, TotalG))+labs(x="Allele Score",y="AUC")+
  geom_jitter(width=0.5,aes(colour = sum.mx1))+theme_minimal()+theme(legend.position='none')
g
g+geom_boxplot(aes(middle=mean(data), color=sum.mx1)) 




#### kinetics ####
setwd("~/Dropbox/Heise/ELISA Antibody/qtls/kinetics")

elisas=aucdata
elisas$day=as.factor(elisas$day)
elisas$day=relevel(elisas$day,"7")

for(i in 1:6)
{
  pheno=abs[i]
  
  iso<-elisas[elisas$isotype == pheno,]
  
  ggplot(iso, aes(day, AUC, colour=day))+geom_point()+facet_wrap(~RIX)+theme_minimal()
  ggsave(filename=paste(pheno,"kinetics_bycross.pdf",sep="_"),width=12,height=10)
}

# rix.list=c(
#   "CC001xCC074","CC004xCC011","CC004xCC012","CC009xCC040",
#   "CC010xCC015","CC011xCC042","CC012xCC041","CC017xCC004",
#   "CC017xCC041","CC018xCC009","CC019xCC002","CC019xCC004",
#   "CC023xCC025","CC024xCC023","CC026xCC034","CC030xCC023",
#   "CC030xCC061","CC032xCC017","CC039xCC020","CC040xCC003",
#   "CC041xCC012","CC042xCC017","CC043xCC037","CC046xCC068",
#   "CC051xCC005","CC061xCC026","CC061xCC039","CC075xCC035")
# 
# iso=iso[which(iso$RIX %in% rix.list),]
# 
# ggplot(iso, aes(day, AUC, colour=day))+geom_point()+facet_wrap(~RIX)+theme_minimal()+
#   labs(x="Day",y="AUC")+theme(legend.position='none')
# ggsave(filename=paste(pheno,"subset_kinetics_bycross.pdf",sep="_"),width=10,height=7)


#### to get into kinetics: strain averages per isotype comparing days ####
library(doBy)

Map.full=NULL
#Map.full aggregates all of the transformed data for each day
#which are different between days - so not directly comparable
#but fine for graphing when separating out days

for(k in 1:4)
{
  day=day.list[[k]]
  dat=dat.list[[k]]
  
  Map.dat=summaryBy(IgG1+IgG2ac+IgG2b+IgG3+IgM+TotalG ~ 
                      RIX + day, data=dat, FUN=mean, na.rm=T)
  Map.full <- rbind(Map.full, Map.dat)
}

colnames(Map.full)[3:8] = gsub(".mean","",colnames(Map.full[3:8]))
Map.full$RIX=as.factor(Map.full$RIX)


# #avg.full is strain averages for full data (all days) w/ no transformations
# 
# avg.full=summaryBy(IgG1+IgG2ac+IgG2b+IgG3+IgM+TotalG ~ 
#                       RIX + day, data=full.dat, FUN=mean, na.rm=T)
# 
# colnames(avg.full)[3:8] = gsub(".mean","",colnames(avg.full[3:8]))
# avg.full$RIX=as.factor(avg.full$RIX)

#graphing that
library(gridExtra)

for(i in 1:6)
{
  phenotype=Map.full[c(1:2,i+2)]
  
  pheno.list=list()
  for (q in 1:3)
  {
    
    phenotype.1=subset(phenotype,phenotype$day==paste(day.list[[q]]))
    colnames(phenotype.1)[3]=paste("D",day.list[[q]],sep="")
    phenotype.1=phenotype.1[c(1,3)]
    phenotype.2=subset(phenotype,phenotype$day==paste(day.list[[q+1]]))
    colnames(phenotype.2)[3]=paste("D",day.list[[q+1]],sep="")
    phenotype.2=phenotype.2[c(1,3)]
    phenotype.3=merge(phenotype.1,phenotype.2,all=T)
    
    g=ggplot(phenotype.3,aes_string(
      x=paste0("D",day.list[[q]],sep=""),
      y=paste0("D",day.list[[q+1]],sep="")
    ))
    gg=g+geom_point(aes(color=RIX),na.rm=T)+theme_bw()+theme(legend.position='none')+
      labs(title=paste(abs[i]))
    
    pheno.list[[q]]=gg
    
  }
  
  pdf(file.path(paste("qtls/kinetics/isotype_kinetics_",abs[i],".pdf",sep="")),width=8,height=3)
  grid.arrange(pheno.list[[1]],pheno.list[[2]],pheno.list[[3]],ncol=3)
  dev.off()
  
}


#### kinetics: two-phenotype isotype splays across days
# Full.trans=rbind(dat.7,dat.10)
# Full.trans=rbind(Full.trans,dat.15)
# Full.trans=rbind(Full.trans,dat.45)

for (k in 1:4)
{
  day=day.list[[k]]
  dat=dat.list[[k]]
  
  pdf(file.path(paste("qtls/d",day,sep=""),paste("pairwise_isotypes_d",day.list[[k]],".pdf",sep="")), width=6, height=6)
  
  pairs(dat[,5:10],lower.panel=panel.smooth,main=paste("Day",day.list[[k]]))
  dev.off()
}


#### comparing total quantities of isotypes ####
library(reshape2)

#day to look at
aucdata=melt(wideauc,id.vars=c("RIX","ID","day","assay_date"),variable.name="isotype",value.name="AUC")

aucdata$day=as.factor(aucdata$day)
auc.plot=subset(aucdata,aucdata$day==10)

ggplot(auc.plot,aes(x=reorder(RIX,AUC), y=AUC,color=isotype))+
  geom_point()+theme_classic()+theme(axis.text.x=element_text(angle=60,vjust=1,hjust=1))

auc.dat=summaryBy(IgG1+IgG2ac+IgG2b+IgG3+IgM+TotalG ~ 
                    RIX + day, data=wideauc, FUN=mean, na.rm=T)
colnames(auc.dat)[3:8] = gsub(".mean","",colnames(wideauc[5:10]))

auc.plot=melt(auc.dat,id.vars=c("RIX","day"),variable.name="isotype",value.name="AUC")
auc.plot=subset(auc.plot,auc.dat$day==10)

ggplot(auc.plot,aes(x=reorder(RIX,AUC), y=AUC,color=isotype))+
  geom_point()+theme_classic()+theme(axis.text.x=element_text(angle=60,vjust=1,hjust=1))

isotype.avgs=summaryBy(IgG1+IgG2ac+IgG2b+IgG3+IgM+TotalG ~ 
                         day, data=auc.dat, FUN=median, na.rm=T)
isotype.avgs=isotype.avgs[1:6] #cuts out totalG
auc.plot=melt(isotype.avgs,id.vars="day",variable.name="isotype",value.name="median.AUC")
ggplot(auc.plot,aes(x=day, y=median.AUC,color=isotype))+
  geom_bar(stat="identity",aes(fill=isotype))+theme_classic()+
  labs(x="Day",y="Median AUC across all strains")

ggsave("isotype_distribution_avg.png",width=6,height=5)