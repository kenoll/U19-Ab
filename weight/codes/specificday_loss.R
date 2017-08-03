setwd("~/Dropbox/Heise/weight_weights")

# raw.weights=read.csv("weights_2016_10.csv")
# 
# raw.weights=raw.weights[,c(1:3,5:7,9:24,40)]
# raw.weights[7:23] = sapply(raw.weights[7:23], as.character)
# raw.weights[7:23] = sapply(raw.weights[7:23], as.numeric)
# 
# raw.weights$D1per=raw.weights$D1/raw.weights$D0*100 
# raw.weights$D2per=raw.weights$D2/raw.weights$D0*100 
# raw.weights$D3per=raw.weights$D3/raw.weights$D0*100 
# raw.weights$D4per=raw.weights$D4/raw.weights$D0*100 
# raw.weights$D5per=raw.weights$D5/raw.weights$D0*100 
# raw.weights$D6per=raw.weights$D6/raw.weights$D0*100 
# raw.weights$D7per=raw.weights$D7/raw.weights$D0*100 
# raw.weights$D8per=raw.weights$D8/raw.weights$D0*100 
# raw.weights$D9per=raw.weights$D9/raw.weights$D0*100
# raw.weights$D10per=raw.weights$D10/raw.weights$D0*100
# raw.weights$D11per=raw.weights$D11/raw.weights$D0*100
# raw.weights$D12per=raw.weights$D12/raw.weights$D0*100
# raw.weights$D13per=raw.weights$D13/raw.weights$D0*100
# raw.weights$D14per=raw.weights$D14/raw.weights$D0*100 
# raw.weights$D15per=raw.weights$D15/raw.weights$D0*100 
# raw.weights$D45per=raw.weights$D45/raw.weights$D0*100 
# 
# weights=raw.weights[c(1:6,24:39)]
# 
# weights=subset(weights,weights$Treatment=="Influenza")
# 
# CC_names=read.csv("~/Dropbox/Heise/CC/cc_names.csv")
# CC_names$Alias=as.character(CC_names$Alias)
# CC_names$CCLine=as.character(CC_names$CCLine)
# 
# colnames(weights)[1]="RIX"
# 
# weights$RIX=as.character(weights$RIX)
# 
# RIX_names <- data.frame(do.call("rbind", strsplit(weights$RIX,"x")))
# 
# RIX_names[,1]=as.character(RIX_names[,1])
# RIX_names[,2]=as.character(RIX_names[,2])
# 
# RIX_CC_1=data.frame(RI_1=RIX_names$X1, CC_1=CC_names[match(RIX_names$X1, CC_names$Alias),1])
# RIX_CC_2=data.frame(RI_2=RIX_names$X2, CC_2=CC_names[match(RIX_names$X2, CC_names$Alias),1])
# 
# RIX_CC=cbind(RIX_CC_1,RIX_CC_2)
# RIX_CC_names=paste(RIX_CC[,2],RIX_CC[,4],sep="x")
# 
# weights[,1]=RIX_CC_names
# 
# write.csv(weights,"weights_pers_2016_10.csv",row.names=F)

weights=read.csv("weights_pers_2016_10.csv")

#day of peak weight weights - returns index of colum with minumum value
weights$Day=as.character(weights$Day)
weights$Day=as.numeric(weights$Day)

####### QTL ANALYSIS ########
library(MASS)

#for D7
weights=weights[c(1:2,4,6,13)]
weights=subset(weights,weights$day>=7)

#set number of phenotype you're looking at
p=1

#set date
run.date="2017_02_20"

#load libraries 
library(ggplot2)

abs=c("D7.wl")
col.list=c("RIX" , "ID" , "cohort" , "day", abs)

colnames(weights)=col.list

#### transformations ####

dat=weights

#remove clear outliers/types (lowest weight <60% starting weight)
# low=which(dat$lowest.weight<60)
# dat=dat[-low,]

for(i in 1:1)
{
  phenotype<-dat[,i+4]
  hist(phenotype, main=abs[i],col="lavender")
}

for(i in 1:1)
{
  phenotype<-dat[,i+4]
  lm.bc=lm(phenotype~RIX,data=dat)
  bc=boxcox(lm.bc)
  print(bc$x[which.max(bc$y)])
}

# Trans.dat.2<-dat[,1:4]
# for(i in 1:p)
# {
#   phenotype<-dat[,i+4]
#   l.phenotype<-sqrt(phenotype)
#   Trans.dat.2<-cbind(Trans.dat.2, l.phenotype)
# }
# colnames(Trans.dat.2)<-col.list


Trans.dat.3<-dat[,1:4]
for(i in 1:1)
{
  phenotype<-dat[,i+4]
  l.phenotype<-(phenotype)^2
  Trans.dat.3<-cbind(Trans.dat.3, l.phenotype)
}
colnames(Trans.dat.3)<-col.list

#lowest day = sqrt
#lowest weight = square
dat=cbind(dat[1:4],Trans.dat.3[5])
dat=dat[col.list]
dat=dat[complete.cases(dat),]

#### phenotypic distribution plots ####
dist.lab=c("day","(percentage of starting weight)^2")

for (i in 1:p)
{
  phenotype=dat[c(1:4,i+4)]
  
  g<-ggplot(phenotype, aes_string(x=abs[i],y=paste0("reorder(RIX,",abs[i],")")))
  g+geom_jitter(aes(color=cohort))+theme_classic()+
    labs(title=paste(abs[i]),y="RIX",x=dist.lab[i])+
    theme(legend.position='none')
  # +xlim(65,150)
  
  ggsave(file.path(paste("qtls/pheno_dist_",abs[i],".jpg",sep="")), width=8, height=8)
}


#### heritability calculations ####
heritability.bounds=NULL
N=3

for(i in 1:p)
{
  a<-anova(lm(dat[,i+4]~dat$RIX))
  interclass.corr<-(a$"Mean Sq"[1] - a$"Mean Sq"[2])/(a$"Mean Sq"[1]+(N-1)*a$"Mean Sq"[2])
  Coef.Genet.Det<-(a$"Mean Sq"[1] - a$"Mean Sq"[2])/(a$"Mean Sq"[1]+((2*N)-1)*a$"Mean Sq"[2])
  h.b=c(interclass.corr,Coef.Genet.Det)
  heritability.bounds <- rbind(heritability.bounds, h.b)
}

rownames(heritability.bounds)=colnames(dat[5:(4+p)])
heritability.bounds=as.data.frame(heritability.bounds)
colnames(heritability.bounds)=c("upper","lower")

write.csv(heritability.bounds,
          file.path(paste("qtls/heritability_estimates","_",abs,".csv",sep="")),
)
heritability.bounds=NULL


#### correlation matrix is meaningless with only two phenotypes ####
#   V<-var(dat[,5:6], na.rm=T)
#   C<-cov2cor(V)
#   
#   jpeg(file.path(paste("qtls/d",day,sep=""),paste("correlation_matrix_d",day.list[[k]],"_",run.date,".jpg",sep="")))
#   heatmap(C, symm=TRUE, col = cm.colors(256))
#   dev.off()

#### mapping ####
#load in required libraries and data
library(doBy)
library(DOQTL)
load("~/Dropbox/Heise/ELISA Antibody/CCRIXb38F.Rdata")
load("~/Dropbox/Heise/ELISA Antibody/MM_snps.Rdata")

model.probs<-model.probs+(1e-20)
K=kinship.probs(model.probs)

## to graph QTL scans (no significance threshold yet)

Map.dat = summaryBy(D7.wl ~ 
                      RIX, data=dat, FUN=mean, na.rm=T)
colnames(Map.dat)[2] = gsub(".mean","",colnames(Map.dat[2]))

Map.dat$Sex="F"
Map.dat$RIX=as.factor(Map.dat$RIX)

row.names(Map.dat)<-Map.dat$RIX
covar = data.frame(sex=as.numeric(Map.dat$Sex == "F"))
rownames(covar)=rownames(Map.dat)


for(i in 1:p)
{
  pheno=abs[i]
  png(file.path(paste("qtls/qtl_scan","_",pheno,"_",run.date,".png",sep="")),width=800,height=400)
  qtl=scanone(pheno=Map.dat, pheno.col=pheno, addcovar=covar, probs=model.probs, K=K, snps=MM_snps)
  plot(qtl,main=paste(pheno))
  #     saveRDS(qtl,file.path(paste("qtls/d",day,sep=""),paste("qtl_scan_d",day,"_",abs[i],"_",run.date,".rds",sep="")))
  save(qtl,file=file.path(paste("qtls/qtl_scan",pheno,".RData",sep="")))
  dev.off()
}


#### with thresholds

for(i in 1:p)
{
  pheno=abs[i]
  png(file.path(paste("qtls/qtl_scan_95","_",pheno,"_",run.date,".png",sep="")),width=800,height=400)
  load(file=file.path(paste("qtls/qtl_scan",pheno,".RData",sep="")))
  
  perms = scanone.perm(pheno = Map.dat, pheno.col = abs[i], probs = model.probs, addcovar = covar, snps= MM_snps, nperm = 100)
  save(perms,file=file.path(paste("qtls/perms100_",abs[i],"_",run.date,".Rdata",sep="")))
  thr = quantile(perms, probs = 0.95)    
  
  plot(qtl,sig.thr=thr,main=paste(pheno))
  dev.off() 
}


thr=quantile(perms, probs = 0.9)
plot(qtl,sig.thr = thr)



# lm.pk=lm(lowest.weight~Treatment+RIX,data=weights)
# summary(lm.pk)
# plot(lm.pk)



## to graph QTL scans w/ Mx1 status as covariate
Map.dat = summaryBy(D7.wl ~ 
                      RIX, data=dat, FUN=mean, na.rm=T)
colnames(Map.dat)[2] = gsub(".mean","",colnames(Map.dat[2]))

#code with Mx1 status - same code as from 2 pheno dist but with score # instead of allele
Map.dat$RIX=as.character(Map.dat$RIX)
RIX_sep <- data.frame(do.call("rbind", strsplit(Map.dat$RIX,"x")))
colnames(RIX_sep)[1:2]=c("dam","sire")
Map.dat=cbind(RIX_sep,Map.dat)
mx1=read.csv("~/Dropbox/Heise/CC/mx1_status.csv")
mx1=mx1[c(1,10)]
colnames(mx1)[1]="dam"
Map.dat=merge(Map.dat,mx1)
colnames(Map.dat)[5]="dam.mx1"
colnames(mx1)[1]="sire"
Map.dat=merge(Map.dat,mx1)
colnames(Map.dat)[6]="sire.mx1"
Map.dat=Map.dat[c(1,6,2,5,3:4)]
Map.dat$mx1.sum=rowSums(Map.dat[c("dam.mx1","sire.mx1")])

Map.dat$Sex="F"
Map.dat$RIX=as.factor(Map.dat$RIX)

#remove cc057 since mx1 haplotype is heterozygous
Map.dat=Map.dat[-which(Map.dat$dam=="CC057" | Map.dat$sire=="CC057"),]

row.names(Map.dat)<-Map.dat$RIX
covar = data.frame(sex=as.numeric(Map.dat$Sex == "F"),mx1=Map.dat$mx1.sum)
rownames(covar)=rownames(Map.dat)
d
for(i in 1:p)
{
  pheno=abs[i]
  png(file.path(paste("qtls/qtl_scan_mx1_covar","_",pheno,"_",run.date,".png",sep="")),width=800,height=400)
  qtl=scanone(pheno=Map.dat, pheno.col=pheno, addcovar=covar, probs=model.probs, K=K, snps=MM_snps)
  plot(qtl,main=paste(pheno))
  #     saveRDS(qtl,file.path(paste("qtls/d",day,sep=""),paste("qtl_scan_d",day,"_",abs[i],"_",run.date,".rds",sep="")))
  save(qtl,file=file.path(paste("qtls/qtl_scan_mx1_covar_",pheno,".RData",sep="")))
  dev.off()
}


for(i in 1:p)
{
  pheno=abs[i]
  png(file.path(paste("qtls/qtl_scan_mx1_covar_95","_",pheno,"_",run.date,".png",sep="")),width=800,height=400)
  load(file=file.path(paste("qtls/qtl_scan_mx1_covar_",pheno,".RData",sep="")))
  
  perms = scanone.perm(pheno = Map.dat, pheno.col = abs[i], probs = model.probs, addcovar = covar, snps= MM_snps, nperm = 100)
  save(perms,file=file.path(paste("qtls/perms100_covar_",abs[i],"_",run.date,".Rdata",sep="")))
  thr = quantile(perms, probs = 0.95)    
  
  plot(qtl,sig.thr=thr,main=paste(pheno))
  dev.off() 
}


### zoom in on allele effects
pheno="d7.wl"
load(file=file.path(paste("qtls/qtl_scan",pheno,".RData",sep="")))

chr=5

coefplot(qtl, chr=chr,legend=F)
coefplot_v2(qtl, chr=chr,main=paste(pheno,"Chromosome",chr,sep=" "),legend=F)

thr = quantile(perms, probs = 0.95) 

#max(qtl$lod$A$lod)
qtl$lod$A[which.max(qtl$lod$A$lod),]
qtl$lod$A[(which.max(qtl$lod$A$lod)-50):(which.max(qtl$lod$A$lod)+20),]


qtl$lod$A[which.max(qtl$lod$A$lod),]
qtl$lod$A[(which.max(qtl$lod$A$lod)-50):(which.max(qtl$lod$A$lod)+20),]


#get confidence interval with 1.5 LOD drop
lod.peak=qtl$lod$A[which.max(qtl$lod$A$lod),"lod"]
lod.min=lod.peak-1.5
sig.range=subset(qtl$lod$A,qtl$lod$A$lod>lod.min & qtl$lod$A$chr==chr)
sig.range[1,]
sig.range[nrow(sig.range),]

#compare to bayesian credible interval
bayesint(qtl,chr=chr,prob=0.95,expandtomarkers=T)

#look at chr5
qtl.chr=subset(qtl$lod$A,qtl$lod$A$chr==5)
qtl.chr[which.max(qtl.chr$lod),]
qtl.chr[(which.max(qtl.chr$lod)-20):(which.max(qtl.chr$lod)+20),]

#to zoom in on qtl plot in region of interest
qtl.chr.sub=subset(qtl.chr,qtl.chr$pos>30 & qtl.chr$pos<55)
# plot(x=qtl.chr.sub$pos,y=qtl.chr.sub$lod)
ggplot(qtl.chr.sub, aes(pos, lod))+geom_point()+theme_minimal()

