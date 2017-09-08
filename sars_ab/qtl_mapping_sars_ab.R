#### setup: load libraries, data format, etc ####
setwd("~/Dropbox/Heise/U19-Ab/sars_ab/")

#load libraries 
library(ggplot2)
library(doBy)
library(DOQTL)

#genome builds
load("~/Dropbox/Heise/CC/cc_refs/CCRIXb38F.Rdata")

#snps info
#load(url("ftp://ftp.jax.org/MUGA/MM_snps.Rdata"))
load("~/Dropbox/Heise/CC/cc_refs/MM_snps.Rdata")
model.probs=model.probs+(1e-20)
K=kinship.probs(model.probs)

#select antigen
anti="S"

#load in transformed data
dataset="lastpos"

dat.7=read.csv(paste0("data_bin/SARS_D7_",anti,"_",dataset,".csv"))
dat.10=read.csv(paste0("data_bin/SARS_D10_",anti,"_",dataset,".csv"))
dat.15=read.csv(paste0("data_bin/SARS_D15_",anti,"_",dataset,".csv"))
dat.29=read.csv(paste0("data_bin/SARS_D29_",anti,"_",dataset,".csv"))
dat.rechal=read.csv(paste0("data_bin/SARS_DRECHAL_",anti,"_",dataset,".csv"))

#set number of ID header columns before data
hc=6

#capture list of column and antibody isotype names as they appear in data
col.list=colnames(dat.10)
abs=colnames(dat.10[(hc+1):length(dat.10)])

#### lists ####
day.list=list("7","10","15","29","rechal")
dat.list=list(dat.7,dat.10,dat.15,dat.29,dat.rechal)

#### phenotypic distribution plots ####

for(k in 1:5)
{
  day=day.list[[k]]
  dat=dat.list[[k]]
  
  for(i in 1:6)
  {
    phenotype=dat[c(1:hc,i+hc)]
    
    g<-ggplot(phenotype, aes_string(x=abs[i],y=paste0("reorder(RIX,",abs[i],")")))
    g+geom_point(na.rm=T)+theme_classic()+
      labs(title=paste("Day",day.list[[k]],abs[i],sep=" "),y="RIX",x=dataset)
    
    ggsave(file.path(paste("phenotypes/d",day,"/pheno_dist",sep=""),paste("pheno_dist_",abs[i],"_",anti,"_",dataset,"_d",day.list[[k]],".jpg",sep="")), width=6, height=9)
  }
}

#### heritability calculations ####
heritability.bounds=NULL
N=3

for(k in 1:5)
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
            file.path(paste("phenotypes/d",day,sep=""),paste("heritability_estimates_",anti,"_",dataset,"_d",day.list[[k]],".csv",sep=""))
  )
  heritability.bounds=NULL
}

#### correlation matrix ####
for(k in 1:5)
{
  day=day.list[[k]]
  dat=dat.list[[k]]
  
  V<-var(dat[(hc+1):length(dat)], na.rm=T)
  C<-cov2cor(V)
  
  jpeg(file.path(paste("phenotypes/d",day,sep=""),paste("correlation_matrix_",anti,"_",dataset,"_d",day.list[[k]],".jpg",sep="")))
  heatmap(C, symm=TRUE, col = cm.colors(256))
  dev.off()
}


#### MAPPING! ####

## to graph QTL scans (no significance threshold yet)
for(k in 1:5)
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
  
  png(file.path(paste("phenotypes/d",day,"/qtl_scans",sep=""),paste("qtl_scan_d",day,"_",anti,"_",dataset,".png",sep=""))
    ,width=1000,height=500)
  par(mfrow = c(2,3))
  for(i in 1:6)
  {
    pheno=abs[i]
    qtl=scanone(pheno=Map.dat, pheno.col=pheno, addcovar=covar, probs=model.probs, K=K, snps=MM_snps)
    plot(qtl,main=paste("Day",day.list[[k]],pheno,sep=" "))
    save(qtl,file=file.path(paste("phenotypes/d",day,"/qtl_scans",sep=""),paste("d",day,"_",abs[i],"_",anti,"_",
                                                                   dataset,"_qtl_obj",".RData",sep="")))
  }
  dev.off()
}

# save PDFs w/ width of 11 and height of 6
# save PNGs w/ width of 1000 and height of 500


#### permutation test to calculate threshold ####

#set number of permutations to run and threshold
#uses previously mapped and saved QTL object
nperm=100
probs=.95

for(k in 1:5)
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
  
  for(i in 1:6)
  {
    pheno=abs[i]
  
    load(file.path(paste("phenotypes/d",day,"/qtl_scans",sep=""),paste("d",day,"_",abs[i],"_",anti,"_",
                                                          dataset,"_qtl_obj",".RData",sep="")))
    
    perms = scanone.perm(pheno = Map.dat, pheno.col = abs[i], probs = model.probs, addcovar = covar, snps= MM_snps, nperm = nperm)
    
    save(perms,
         file=file.path(paste("phenotypes/d",day,"/qtl_scans",sep=""),
                   paste("d",day,"_",abs[i],"_",anti,"_",dataset,"_perms_",nperm,".Rdata",sep="")))
    
    thr = quantile(perms, probs = probs)    
    
    plot(qtl,sig.thr=thr,main=paste("Day",day,anti,pheno,dataset,sep=" "))
  }
  dev.off()
}

### zoom in on chromosomes with significant/suggestive QTL to look at range and allele effects ###
day=7
pheno="TotalG"
chr=14

#load QTL object
load(file=file.path(paste("phenotypes/d",day,"/qtl_scans",sep=""),paste("d",day,"_",pheno,"_",anti,"_",
                                                      dataset,"_qtl_obj",".RData",sep="")))
coefplot(qtl, chr=chr,legend=F)

### to get regions of genome that fall under significance threshold
#load permutations
load(file=file.path(paste("phenotypes/d",day,"/qtl_scans",sep=""),
               paste("d",day,"_",pheno,"_",anti,"_",dataset,"_perms_",nperm,".Rdata",sep="")))
thr = quantile(perms, probs = 0.95) 

#get marker with maximum LOD score and look at surrounding interval
qtl$lod$A[which.max(qtl$lod$A$lod),]
# qtl$lod$A[(which.max(qtl$lod$A$lod)-10):(which.max(qtl$lod$A$lod)+10),]

#get confidence interval with 1.5 LOD drop
lod.peak=qtl$lod$A[which.max(qtl$lod$A$lod),"lod"]
lod.min=lod.peak-1.5
sig.range=subset(qtl$lod$A,qtl$lod$A$lod>lod.min & qtl$lod$A$chr==chr)
sig.range[1,] #lower end of confidence interval
sig.range[nrow(sig.range),] #upper end of confidence interval

#compare to bayesian credible interval
bayesint(qtl,chr=chr,prob=0.95,expandtomarkers=T)


#### more fun ####

#get haplotype probs for RIXs at the peak - manually 
model.probs=model.probs-(1e-20)
hap.probs<-model.probs[,,qtl$lod$A[which.max(qtl$lod$A$lod),"marker"]]
# write.csv(hap.probs,file.path(paste("phenotypes/d",day,sep=""),
#               paste("d",day,"_",pheno,"_",anti,"_",dataset,"_chr",chr,"_hap_probs",".csv",sep="")))

#to look at a peak other than the very max, input chromosome of interest
chr2=8
qtl.chr=subset(qtl$lod$A,qtl$lod$A$chr==chr2)
qtl.chr[which.max(qtl.chr$lod),]

#to zoom in on qtl plot in region of interest (e.g. if there are weird double peaks)
qtl.chr.sub=subset(qtl.chr,qtl.chr$pos>100 & qtl.chr$pos<120)
ggplot(qtl.chr.sub, aes(pos, lod))+geom_point()+theme_minimal()
