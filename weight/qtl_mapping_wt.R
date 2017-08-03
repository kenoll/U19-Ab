#### Setup ####
setwd("~/Dropbox/Heise/U19-Ab/weight/")

#load libraries 
library(ggplot2)
library(dplyr)
library(DOQTL)
library(doBy)

#load in transformed data (see "Antibody data for mapping pipeline.R")
dat=read.csv("weights_pers_2016_10.csv")
  dat=dat[c(1:6,13)]
  dat=dat[complete.cases(dat),]
  colnames(dat)[7]="D7"
dat2=read.csv("peak_weight_loss.csv")
  dat2=dat2[-7]
dat=merge(dat,dat2)
  rm(dat2)
    
save.as="wl"

#set number of ID columns in phenotype sheet and QTL chromosome
hc=6
abs=colnames(dat)[(hc+1):length(dat)]
col.list=colnames(dat)

#load in required data
load("~/Dropbox/Heise/U19-Ab/cc_refs/CCRIXb38F.Rdata") #RIX model probs
# load("~/Dropbox/Heise/U19-Ab/cc_refs/CCModelProbs.Rdata") #RI model probs
# load(url("ftp://ftp.jax.org/MUGA/MM_snps.Rdata"))
load("~/Dropbox/Heise/U19-Ab/cc_refs/MM_snps.Rdata")

model.probs<-model.probs+(1e-20)
K=kinship.probs(model.probs)

## to graph QTL scans (w/ or w/o Mx1 status as covariate)
Map.dat=dat[c(1,(hc+1):(hc+length(abs)))]
Map.dat=summaryBy(. ~ RIX, data=dat, FUN=mean, na.rm=T)
colnames(Map.dat)[2:(1+length(abs))] = gsub(".mean","",colnames(Map.dat[2:(1+length(abs))]))
 
#code with Mx1 status - see add.stat.mx1.R for function code
Map.dat=add.stat.mx1(Map.dat)
      
#create covariate information
row.names(Map.dat)<-Map.dat$RIX
Map.dat$Sex="F"
covar.sex = data.frame(sex=as.numeric(Map.dat$Sex == "F"))
covar.mx1 = data.frame(sex=as.numeric(Map.dat$Sex == "F"),mx1=Map.dat$sum.mx1)
rownames(covar.sex)=rownames(Map.dat)
rownames(covar.mx1)=rownames(Map.dat)
      


#### Basic QTL plot ####

  #choose whether you want mx1 as a covariate or not
  covar.stat="mx1"
  # covar.stat=F
  
  if (covar.stat=="mx1") {
    covar=covar.mx1;covar.save="_mx1_covar"} else {
      covar=covar.sex;covar.save=""}
  
  fp=paste0("qtl_scans/qtl_scan_",save.as,covar.save)
  
  #generate QTL object and save image of scan
    for(i in 1:length(abs))
    {
      pheno=abs[i]
      png(file.path(paste0(fp,"_",pheno,".png")),width=800,height=400)
      qtl=scanone(pheno=Map.dat, pheno.col=pheno, addcovar=covar, probs=model.probs, K=K, snps=MM_snps)
      plot(qtl,main=paste(pheno))
      save(qtl,file=file.path(paste0(fp,"_",pheno,".RData")))
      dev.off()
    }
  
  
#### Significance Testing ####

  #set number of permutations to run and threshold
  nperm=100
  probs=.95

  #loads in QTL object from above (rather than generating a new one)
  #does permutations, saves permutations, calculates significance threshold
  #creates new plot with significance threshold at the specified confidence interval
   for(i in 1:length(abs))
   {
      pheno=abs[i]
      png(file.path(paste0(fp,"_",pheno,"_thr",probs*100,".png")),width=800,height=400)
      load(file.path(paste0(fp,"_",pheno,".RData")))
      perms = scanone.perm(pheno = Map.dat, pheno.col = pheno, probs = model.probs, addcovar = covar, snps= MM_snps, nperm = nperm)
      save(perms,
            file=file.path(paste0(fp,"_",pheno,"_perms_",nperm,".Rdata")))
      thr = quantile(perms, probs = probs)
      plot(qtl,sig.thr=thr,main=paste(pheno))
   }
      dev.off()

      
#### Analyze peak ####
  
  # select phenotype of interest
  pheno="D7"
      
  # load in QTL objects
  load(file.path(paste0(fp,"_",pheno,".RData")))

  # confirm QTL and select chromosome of interest
  plot(qtl)
  chr=1
  
  # generate allele effects plot
  coefplot(qtl, chr=chr,legend=F)
  #coefplot_v2(qtl, chr=chr,main=paste(pheno,"Chromosome",chr,sep=" "),legend=F)

  ### get coordinates for regions under significance thresold
  load(file=file.path(paste0(fp,"_",pheno,"_perms_",nperm,".Rdata")))
  thr = quantile(perms, probs = probs)

  #get marker with maximum LOD score and look at surrounding interval
  qtl$lod$A[which.max(qtl$lod$A$lod),]
  # qtl$lod$A[(which.max(qtl$lod$A$lod)-10):(which.max(qtl$lod$A$lod)+10),]

  #get confidence interval with 1.5 LOD drop
  lod.peak=qtl$lod$A[which.max(qtl$lod$A$lod),"lod"]
  lod.min=lod.peak-1.5
  sig.range=subset(qtl$lod$A,qtl$lod$A$lod>lod.min & qtl$lod$A$chr==chr)
  sig.range[1,] ;  sig.range[nrow(sig.range),]

  # compare to bayesian credible interval
  bayesint(qtl,chr=chr,prob=0.95,expandtomarkers=T)

  # get haplotype probs for RIXs at the peak - manually
  # model.probs[,,qtl$lod$A[which.max(qtl$lod$A$lod),"marker"]])

  
#### Examining Weeds ####  
  
  #to get markers for a peak other than the very max, input chromosome of interest
  chr2 = 4
  qtl.chr=subset(qtl$lod$A,qtl$lod$A$chr==chr2)
  qtl.chr[which.max(qtl.chr$lod),]
  # qtl.chr[(which.max(qtl.chr$lod)-10):(which.max(qtl.chr$lod)+10),]
  
  #to zoom in on qtl plot in region of interest, input regions
    qtl.chr %>% subset(qtl.chr$pos>120 & qtl.chr$pos<150) %>%
      ggplot(aes(pos, lod))+geom_point()+theme_minimal()
  