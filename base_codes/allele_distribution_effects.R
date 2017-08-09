#2016_11_23 updated to used average scores that encompass recombinations within locus

setwd("~/Dropbox/Heise/ELISA Antibody/")
chrom=read.csv("qtls/d10/snps/chr5_68-75_scores.csv")

##### strain averages #####
dat=dat.10

Map.dat=summaryBy(IgG1+IgG2ac+IgG2b+IgG3+IgM+TotalG ~ 
                    RIX + day, data=dat, FUN=mean, na.rm=T)
colnames(Map.dat)[3:8] = gsub(".mean","",colnames(Map.dat[3:8]))
Map.dat$RIX=as.character(Map.dat$RIX)


RIX_sep <- data.frame(do.call("rbind", strsplit(Map.dat$RIX,"x")))
colnames(RIX_sep)[1:2]=c("dam","sire")

Map.dat=cbind(RIX_sep,Map.dat)

chrom=chrom[c(1,length(chrom))]

colnames(chrom)[1]="dam"
Map.dat=merge(Map.dat,chrom)
colnames(Map.dat)[11]="dam.chrom"

colnames(chrom)[1]="sire"
Map.dat=merge(Map.dat,chrom)
colnames(Map.dat)[12]="sire.chrom"

Map.dat=Map.dat[c(1,12,2,11,3:10)]

Map.dat$additive.status=rowSums(Map.dat[c("dam.chrom","sire.chrom")])

Map.dat$additive.status=as.numeric(Map.dat$additive.status)

hist(Map.dat$additive.status,main=paste("Distribution of Founder Alleles at Chr5 QTL"),
     xlab=paste("Allele Score at Chr5 QTL"),col="moccasin")

# png("qtls/d10/chr5_distribution.png",width=500,height=400)
# 
# hist(Map.dat$additive.status,main=paste("Distribution of Founder Alleles at Chr5 QTL"),
#      xlab=paste("Allele Score at Chr5 QTL"),col="moccasin")
# 
# dev.off()

Map.dat$additive.status=as.character(Map.dat$additive.status)
Map.dat$additive.status=as.factor(Map.dat$additive.status)

ggplot(Map.dat, aes(additive.status, IgM))+labs(x="Allele Score at Chr5 QTL",y="Log10 IgM AUC")+
  geom_jitter(width=0.5,aes(colour = additive.status))+theme_minimal()+theme(legend.position='none')

ggsave("qtls/d10/chrom_allele_igm_strainavg.png")



#####mapping with an allele effect as a covariate
# library(doBy)
# library(DOQTL)
# load("~/Dropbox/Heise/ELISA Antibody/R codes/CCRIXb38F.Rdata")
# load(url("ftp://ftp.jax.org/MUGA/MM_snps.Rdata"))
# 
# model.probs<-model.probs+(1e-20)
# K=kinship.probs(model.probs)

  Map.dat$Sex="F"
  Map.dat$RIX=as.factor(Map.dat$RIX)
  Map.dat$additive.status=as.character(Map.dat$additive.status)
  Map.dat$additive.status=as.numeric(Map.dat$additive.status)
  
  row.names(Map.dat)<-Map.dat$RIX
  covar = data.frame(sex=as.numeric(Map.dat$Sex == "F"),chrom=Map.dat$additive.status)
  rownames(covar)=rownames(Map.dat)
  
  png(file.path(paste("qtls/d",day,sep=""),paste("qtl_IgM_chromcovariate_2_d",day.list[[k]],".png",sep="")),width=1000,height=500)
  par(mfrow = c(2,3))
  for(i in 1:6)
  {
    pheno=abs[i]
    qtl=scanone(pheno=Map.dat, pheno.col=pheno, addcovar=covar, probs=model.probs, K=K, snps=MM_snps)
    plot(qtl,main=paste("Day",day.list[[k]],pheno,"w/ chrom as covariate",sep=" "))
    #     saveRDS(qtl,file.path(paste("qtls/d",day,sep=""),paste("qtl_scan_d",day,"_",abs[i],".rds",sep="")))
    #     save(qtl,file=file.path(paste("qtls/d",day,sep=""),paste("qtl_scan_d",day,"_",abs[i],".RData",sep="")))
  }
  dev.off()
}

qtl.chr=subset(qtl$lod$A,qtl$lod$A$chr==5)
qtl.chr.sub=subset(qtl.chr,qtl.chr$pos>63 & qtl.chr$pos<75)
ggplot(qtl.chr.sub, aes(pos, lod))+geom_point()+theme_minimal()



#### individual mice ####
setwd("~/Dropbox/Heise/ELISA Antibody/")
chrom=read.csv("qtls/d10/chr5_status_curated.csv")


dat=dat.10

dat$RIX=as.character(dat$RIX)
RIX_sep <- data.frame(do.call("rbind", strsplit(dat$RIX,"x")))
colnames(RIX_sep)[1:2]=c("dam","sire")

dat=cbind(RIX_sep,dat)

chrom=chrom[c(1,8)]

colnames(chrom)[1]="dam"
dat=merge(dat,chrom)
colnames(dat)[(length(dat))]="dam.chrom"

colnames(chrom)[1]="sire"
dat=merge(dat,chrom)
colnames(dat)[(length(dat))]="sire.chrom"

dat=dat[c(1,14,2,13,3:11)]

dat$additive.status=rowSums(dat[c("dam.chrom","sire.chrom")])

dat$additive.status=as.numeric(dat$additive.status)

hist(dat$additive.status,main=paste("Distribution of Founder Alleles at Chr5 QTL"),
     xlab=paste("Allele Score at Chr5 QTL"),col="moccasin")

# png("qtls/d10/chr5_distribution.png",width=500,height=400)
# 
# hist(dat$additive.status,main=paste("Distribution of Founder Alleles at Chr5 QTL"),
#      xlab=paste("Allele Score at Chr5 QTL"),col="moccasin")
# 
# dev.off()

dat$additive.status=as.character(dat$additive.status)
dat$additive.status=as.factor(dat$additive.status)

ggplot(dat, aes(additive.status, IgG3))+labs(x="Allele Score at Chr5 QTL",y="Log10 IgM AUC")+
  geom_jitter(width=0.5,aes(colour = additive.status))+theme_minimal()+theme(legend.position='none')

ggsave("qtls/d10/chrom_allele_igm.png")

