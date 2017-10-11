#### setup: load libraries, data format, etc ####
setwd("~/Dropbox/Heise/U19-Ab/antibody/")

#load libraries 
library(ggplot2)
library(DOQTL)
library(doBy)
library(gridExtra)
library(reshape2)

#load in transformed data (see "Antibody data for mapping pipeline.R")

data.type="auc"

dat.7=read.csv(paste("data_bin/d7/dat.7.",data.type,".csv",sep=""))
dat.10=read.csv(paste("data_bin/d10/dat.10.",data.type,".csv",sep=""))
dat.15=read.csv(paste("data_bin/d15/dat.15.",data.type,".csv",sep=""))
dat.45=read.csv(paste("data_bin/d45/dat.45.",data.type,".csv",sep=""))

#set number of header columns before data starts
hc=10

#set list of antibody isotypes in the column order they appear in
abs=c("IgG1" ,  "IgG2ac" , "IgG2b" , "IgG3" ,  "IgM"  ,  "TotalG")
col.list=colnames(dat.10)


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
    phenotype=dat[c(1:hc,i+hc)]
    phenotype=phenotype[complete.cases(phenotype),]  
    
    g<-ggplot(phenotype, aes_string(x=abs[i],y=paste0("reorder(RIX,",abs[i],")")))
    g+geom_point(aes(colour=assay_date))+theme_classic()+
      labs(title=paste("Day",day.list[[k]],abs[i],sep=" "),y="RIX",x=data.type)
    
    ggsave(file.path(paste("phenotypes/d",day,"/pheno_dist",sep=""),paste("pheno_dist_",abs[i],"_",data.type,"_d",day.list[[k]],".jpg",sep="")), width=6, height=9)
  }
}

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
    a<-anova(lm(dat[,i+hc]~dat$RIX))
    interclass.corr<-(a$"Mean Sq"[1] - a$"Mean Sq"[2])/(a$"Mean Sq"[1]+(N-1)*a$"Mean Sq"[2])
    Coef.Genet.Det<-(a$"Mean Sq"[1] - a$"Mean Sq"[2])/(a$"Mean Sq"[1]+((2*N)-1)*a$"Mean Sq"[2])
    h.b=c(interclass.corr,Coef.Genet.Det)
    heritability.bounds <- rbind(heritability.bounds, h.b)
  }
  
  rownames(heritability.bounds)=abs
  heritability.bounds=as.data.frame(heritability.bounds)
  heritability.bounds[,3]=rep.int(day.list[[k]],nrow(heritability.bounds))
  colnames(heritability.bounds)=c("upper","lower","day")
  
  write.csv(heritability.bounds,
            file.path(paste0("phenotypes/d",day),
                      paste0("heritability_estimates_",data.type,"_d",day.list[[k]],".csv")))
  heritability.bounds=NULL
}

#### correlation matrix ####
for(k in 1:4)
{
  day=day.list[[k]]
  dat=dat.list[[k]]
  
  V<-var(dat[(hc+1):length(dat)], na.rm=T)
  C<-cov2cor(V)
  
  jpeg(file.path(paste0("phenotypes/d",day),paste("correlation_matrix_",data.type,"_d",day.list[[k]],".jpg",sep="")))
  heatmap(C, symm=TRUE, col = cm.colors(256))
  dev.off()
}

          # #### two-phenotype distributions ####
          # library(doBy)
          # dat=dat.10
          # 
          # Map.dat=summaryBy(IgG1+IgG2ac+IgG2b+IgG3+IgM+TotalG ~ 
          #                     RIX + day, data=dat, FUN=mean, na.rm=T)
          # colnames(Map.dat)[3:8] = gsub(".mean","",colnames(Map.dat[3:8]))
          # Map.dat$RIX=as.character(Map.dat$RIX)
          # 
          # #code with Mx1 status
          # 
          # RIX_sep <- data.frame(do.call("rbind", strsplit(Map.dat$RIX,"x")))
          # colnames(RIX_sep)[1:2]=c("dam","sire")
          # 
          # Map.dat=cbind(RIX_sep,Map.dat)
          # 
          # mx1=read.csv("~/Dropbox/Heise/CC/mx1_status.csv")
          # mx1=mx1[c(1,8)]
          # 
          # colnames(mx1)[1]="dam"
          # Map.dat=merge(Map.dat,mx1)
          # colnames(Map.dat)[11]="dam.mx1"
          # 
          # colnames(mx1)[1]="sire"
          # Map.dat=merge(Map.dat,mx1)
          # colnames(Map.dat)[12]="sire.mx1"
          # 
          # Map.dat=Map.dat[c(1,12,2,11,3:10)]
          # 
          # #graphs
          # Map.dat$RIX=as.factor(Map.dat$RIX)
          # 
          # phenotype=Map.dat[c(1:6,12,11)]
          # phenotype=phenotype[complete.cases(phenotype),]  
          # 
          # cbPalette <- c("royalblue","orchid","seagreen2")
          # 
          # g=ggplot(phenotype,aes(x=TotalG,y=IgM))
          # g+geom_point(size=4,aes(color=dam.mx1))+
          #   scale_color_manual(values=cbPalette)+
          #   geom_point(aes(color=sire.mx1))+
          #   theme_classic()+geom_text(aes(label=RIX),size=3,vjust=-1)+
          #   theme(legend.title=element_blank())
          # 
          # 
          # # rix.list=c("CC037xCC046","CC063xCC002","CC051xCC009","CC034xCC016",
          # #            "CC016xCC038","CC016xCC061","CC060xCC006")
          # 
          # # rix.list=("CC074xCC062")
          # 
          # 
          # rix.list=c("CC003xCC062")
          # 
          # rix.subset=Map.dat[which(Map.dat$RIX %in% rix.list),]
          # 
          # g=ggplot(phenotype,aes(x=TotalG,y=IgM))
          # g+geom_point(size=4,aes(color=dam.mx1))+
          #   scale_color_manual(values=cbPalette)+
          #   geom_point(aes(color=sire.mx1))+
          #   theme_classic()+geom_text(aes(label=RIX),size=3,vjust=-1)+
          #   theme(legend.title=element_blank())+
          #   geom_point(data=rix.subset,aes(x=TotalG,y=IgM),colour="red", size=5)
          # 
          # 
          # # ggsave("qtls/d45/totalg_v_igm_d7_highlighted.jpg",width=10,height=6)


#### mapping ####
#load in required data
load("~/Dropbox/Heise/CC/cc_refs/CCRIXb38F.Rdata")
# load(url("ftp://ftp.jax.org/MUGA/MM_snps.Rdata"))
load("~/Dropbox/Heise/CC/cc_refs/MM_snps.Rdata")

model.probs<-model.probs+(1e-20)
K=kinship.probs(model.probs)

## to graph QTL scans (no significance threshold yet)

for(k in 1:4)
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
  
  png(file.path(paste0("phenotypes/d",day,"/qtl_scans/qtl_scan_d",day,"_",data.type,".png"))
    ,width=1000,height=500)
  par(mfrow = c(2,3))
  for(i in 1:6)
  {
    pheno=abs[i]
    qtl=scanone(pheno=Map.dat, pheno.col=pheno, addcovar=covar, probs=model.probs, K=K, snps=MM_snps)
    plot(qtl,main=paste("Day",day.list[[k]],pheno,sep=" "))
    save(qtl,file=file.path(paste0("phenotypes/d",day,"/qtl_scans/d",day,"_",abs[i],"_",
                                   data.type,"_qtl_obj",".RData")))
  }
  dev.off()
}

# save PDFs w/ width of 11 and height of 6
# save PNGs w/ width of 1000 and height of 500

#### permutation test to calculate threshold ####
#set number of permutations to run and threshold
nperm=100

for(k in 1:4)
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
  

  png(file.path(paste0("phenotypes/d",day,"/qtl_scans/qtl_scan_d",day,"_",data.type,
                           "_perms_",nperm,".png")) ,width=1000,height=500)
  par(mfrow = c(2,3))
  for(i in 1:6)
  {
    pheno=abs[i]
  
    load(file.path(paste0("phenotypes/d",day,"/qtl_scans/d",day,"_",abs[i],"_",
                     data.type,"_qtl_obj",".RData")))
    
    perms = scanone.perm(pheno = Map.dat, pheno.col = abs[i], probs = model.probs, addcovar = covar, snps= MM_snps, nperm = nperm)
    
    save(perms,
         file=file.path(paste0("phenotypes/d",day,"/qtl_scans/d",day,"_",abs[i],"_",
                  data.type,"_perms_",nperm,".Rdata",sep="")))
    
    thr = quantile(perms, probs = 0.95) 
    
    thr2 = quantile(perms, probs = 0.9)  
    
    plot(qtl,sig.thr=thr,main=paste("Day",day,pheno,data.type,sep=" "))
    
    abline(h=thr2,col="blue")    
  }
  dev.off()
}

### zoom in on chromosomes with significant/suggestive QTL to look at range and allele effects ###
day=7
pheno="IgG2ac"
load(file=file.path(paste0("phenotypes/d",day,"/qtl_scans/d",day,"_",pheno,"_",
                                                      data.type,"_qtl_obj",".RData")))

chromo=17

# plot(qtl, main=paste("Day",day,pheno,sep=" "))
coefplot(qtl, chr=chromo,main=paste("Day",day,pheno,"Chromosome",chromo,sep=" "),legend=F)
coefplot_v2(qtl, chr=chromo,main=paste("Day",day,pheno,"Chromosome",chromo,sep=" "),legend=F)

### to get regions ofgenome that fall under significance threshold
load(file=file.path(paste0("phenotypes/d",day,"/qtl_scans/d",day,"_",pheno,"_",
                           data.type,"_perms_",nperm,".Rdata")))

thr = quantile(perms, probs = 0.95) 

#get marker with maximum LOD score and look at surrounding interval
qtl$lod$A[which.max(qtl$lod$A$lod),]
# qtl$lod$A[(which.max(qtl$lod$A$lod)-10):(which.max(qtl$lod$A$lod)+10),]

#get confidence interval with 1.5 LOD drop
lod.peak=qtl$lod$A[which.max(qtl$lod$A$lod),"lod"]
lod.min=lod.peak-1.5
sig.range=subset(qtl$lod$A,qtl$lod$A$lod>lod.min & qtl$lod$A$chr==chromo)
sig.range[c(1,nrow(sig.range)),]

#compare to bayesian credible interval
bayesint(qtl,chr=chromo,prob=0.95,expandtomarkers=T)


#qqnorm plots
x=qtl$lod$A$lod
  sub.x=qtl$lod$A
    sub.x=subset(sub.x,sub.x$chr==chromo)
    sub.x=(sub.x$lod)
  sub.not=qtl$lod$A
    sub.not=subset(sub.not,sub.not$chr!=chromo)
    sub.not=(sub.not$lod)


png(file.path(paste0("phenotypes/d",day,"_",pheno,"_",data.type,"_qqnorm.png")))
  qqnorm(x,col="darkgray",main=paste0("Day",day,pheno),xlim=c(-4,4),ylim=c(0,7));qqline(x)
  dev.off()

png(file.path(paste0("phenotypes/d",day,"_",pheno,"_",data.type,"_qqnorm_chr",chromo,".png")))
  qqnorm(sub.x,col="darkgray",main=paste0("Day",day,pheno,"- Chromosome",chromo),xlim=c(-4,4),ylim=c(0,7));qqline(x)
  dev.off()
  
png(file.path(paste0("phenotypes/d",day,"_",pheno,"_",data.type,"_qqnorm_notchr",chromo,".png")))
  qqnorm(sub.not,col="darkgray",main=paste0("Day",day,pheno,"- All except Chromosome",chromo),xlim=c(-4,4),ylim=c(0,7));qqline(x)
  dev.off()




#get haplotype probs for RIXs at the peak - manually 
probs<-model.probs[,,qtl$lod$A[which.max(qtl$lod$A$lod),"marker"]]
 # write.csv(probs,file.path(paste("phenotypes/d",day,sep=""),
 #          paste("d",day,"_",pheno,"_",data.type,"_chr",chromo,"_probs",".csv",sep="")))


#to look at a peak other than the very max, input chromosome of interest
qtl.chr=subset(qtl$lod$A,qtl$lod$A$chr==5)
qtl.chr[which.max(qtl.chr$lod),]
# qtl.chr[(which.max(qtl.chr$lod)-10):(which.max(qtl.chr$lod)+10),]

#to zoom in on qtl plot in region of interest
qtl.chr.sub=subset(qtl.chr,qtl.chr$pos>63 & qtl.chr$pos<75)
# plot(x=qtl.chr.sub$pos,y=qtl.chr.sub$lod)
ggplot(qtl.chr.sub, aes(pos, lod))+geom_point()+theme_minimal()










#####where i am


## to graph QTL scans w/ Mx1 status as covariate
dat=dat.10
Map.dat=summaryBy(IgG1+IgG2ac+IgG2b+IgG3+IgM+TotalG ~ 
                    RIX + day, data=dat, FUN=mean, na.rm=T)
colnames(Map.dat)[3:8] = gsub(".mean","",colnames(Map.dat[3:8]))
Map.dat$RIX=as.character(Map.dat$RIX)

#code with Mx1 status - same code as from 2 pheno dist but with score # instead of allele
Map.dat$RIX=as.character(Map.dat$RIX)
RIX_sep <- data.frame(do.call("rbind", strsplit(Map.dat$RIX,"x")))
colnames(RIX_sep)[1:2]=c("dam","sire")
Map.dat=cbind(RIX_sep,Map.dat)
mx1=read.csv("~/Dropbox/Heise/CC/mx1_status.csv")
mx1=mx1[c(1,10)]
colnames(mx1)[1]="dam"
Map.dat=merge(Map.dat,mx1)
colnames(Map.dat)[11]="dam.mx1"
colnames(mx1)[1]="sire"
Map.dat=merge(Map.dat,mx1)
colnames(Map.dat)[12]="sire.mx1"
# Map.dat=Map.dat[c(1,12,2,11,3:10)]
Map.dat$mx1.sum=rowSums(Map.dat[c("dam.mx1","sire.mx1")])

Map.dat$Sex="F"
Map.dat$RIX=as.factor(Map.dat$RIX)

#remove cc057 since mx1 haplotype is heterozygous
Map.dat=Map.dat[-which(Map.dat$dam=="CC057" | Map.dat$sire=="CC057"),]

row.names(Map.dat)<-Map.dat$RIX
covar = data.frame(sex=as.numeric(Map.dat$Sex == "F"),mx1=Map.dat$mx1.sum)
rownames(covar)=rownames(Map.dat)


# for(i in 1:6)
{
  fp=paste("qtls/qtl_scan_mx1_covar_d",day.list[k],"_",pheno,"_",run.date,sep="")
  pheno=abs[i]
  png(file.path(paste(fp,".png",sep="")),width=800,height=400)
  qtl=scanone(pheno=Map.dat, pheno.col=pheno, addcovar=covar, probs=model.probs, K=K, snps=MM_snps)
  plot(qtl,main=paste(pheno))
  #     saveRDS(qtl,file.path(paste("qtls/d",day,sep=""),paste("qtl_scan_d",day,"_",abs[i],"_",run.date,".rds",sep="")))
  save(qtl,file=file.path(paste(fp,".RData",sep="")))
  dev.off()
}










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
#day to look at
wideauc=wideauc[c(1:3,7:13)]
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
colnames(isotype.avgs)[2:6] = gsub(".median","",colnames(isotype.avgs[2:6]))

isotype.avgs$day=as.factor(isotype.avgs$day)
isotype.avgs$day=relevel(isotype.avgs$day,"7")

auc.plot=melt(isotype.avgs,id.vars="day",variable.name="isotype",value.name="median.AUC")

ggplot(auc.plot,aes(x=day, y=median.AUC,color=isotype))+
  geom_bar(stat="identity",aes(fill=isotype))+theme_minimal()+
  labs(x="Day",y="Median AUC across all strains")+theme(legend.position='bottom')+theme(legend.title=element_blank())

ggsave("isotype_distribution_avg.png",width=6,height=5)
