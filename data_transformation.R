#### setup: load libraries, data format, etc ####
setwd("~/Dropbox/Heise/ELISA Antibody/")

#load libraries 
library(ggplot2)
library(MASS)

##### AUC #####

#load in data to map in wide format
load("wideauc.Rdata")
wideauc$RIX=as.factor(wideauc$RIX)
full.dat<-wideauc

#set list of anitbody isotypes in the column order they appear in
abs=c("IgG1" ,  "IgG2ac" , "IgG2b" , "IgG3" ,  "IgM"  ,  "TotalG")

col.list=c("RIX" , "ID" , "day" , "assay_date", "IgG1" ,"IgG2ac","IgG2b" ,    
           "IgG3","IgM","TotalG")


###### visualize data and transform as necessary (individually by timepoint) #####
## plot data distribution and transform as necessary

#### day 7 ####

dat=subset(full.dat,full.dat$day==7)

# graph to view data distribution #
par(mfrow = c(2,3))

for(i in 1:6)
{
  phenotype<-dat[,i+4]
  hist(phenotype, main=abs[i])
}

#some AUC values are every so slightly below 0.. box cox needs all positive values
#change the values to be positive non-zero numbers so box-cox will run

for(i in 1:6)
{
  phenotype<-dat[,i+4]
  phenotype[phenotype<0] <- 0.00000000001
  lm.bc=lm(phenotype~RIX,data=dat)
  bc=boxcox(lm.bc)
  print(bc$x[which.max(bc$y)])
}

#see what they look like log transformed - more normal?
Trans.dat<-dat[,1:4]

for(i in 1:6)
{
  phenotype<-dat[,i+4]
  l.phenotype<-log(phenotype, base=10)
  Trans.dat<-cbind(Trans.dat, l.phenotype)
}
colnames(Trans.dat)<-c("RIX","ID","day","assay_date","IgG1","IgG2ac","IgG2b","IgG3","IgM","TotalG")

for(i in 1:6)
{
  phenotype<-Trans.dat[,i+4]
  hist(phenotype, main=abs[i])
}


#see what they look like sqrt transformed - more normal?
Trans.dat.2<-dat[,1:4]

for(i in 1:6)
{
  phenotype<-dat[,i+4]
  l.phenotype<-sqrt(phenotype)
  Trans.dat.2<-cbind(Trans.dat.2, l.phenotype)
}
colnames(Trans.dat.2)<-c("RIX","ID","day","assay_date","IgG1","IgG2ac","IgG2b","IgG3","IgM","TotalG")

for(i in 1:6)
{
  phenotype<-Trans.dat.2[,i+4]
  hist(phenotype, main=abs[i])
}


##D7 Transformations
#log transform all
dat.7=cbind(Trans.dat[c(1:10)])
dat.7=dat.7[col.list]

#visualize final data
for(i in 1:6)
{
  phenotype<-dat.7[,i+4]
  hist(phenotype, main=abs[i],col="lavender")
}

write.csv(dat.7,file=paste("~/Dropbox/Heise/ELISA Antibody/qtls/d7/dat.7.auc.csv",sep=""),row.names=F)


#### day 10 ####
dat=subset(full.dat,full.dat$day==10)

#graph to view data distribution
par(mfrow = c(2,3))

for(i in 1:6)
{
  phenotype<-dat[,i+4]
  hist(phenotype, main=abs[i])
}

for(i in 1:6)
{
  phenotype<-dat[,i+4]
  phenotype[phenotype<0] <- 0.00000000001
  lm.bc=lm(phenotype~RIX,data=dat)
  bc=boxcox(lm.bc)
  print(bc$x[which.max(bc$y)])
}


#see what they look like log transformed - more normal?
Trans.dat<-dat[,1:4]

for(i in 1:6)
{
  phenotype<-dat[,i+4]
  l.phenotype<-log(phenotype, base=10)
  Trans.dat<-cbind(Trans.dat, l.phenotype)
}
colnames(Trans.dat)<-c("RIX","ID","day","assay_date","IgG1","IgG2ac","IgG2b","IgG3","IgM","TotalG")


for(i in 1:6)
{
  phenotype<-Trans.dat[,i+4]
  hist(phenotype, main=abs[i])
}


#see what they look like sqrt transformed - more normal?
Trans.dat.2<-dat[,1:4]

for(i in 1:6)
{
  phenotype<-dat[,i+4]
  l.phenotype<-sqrt(phenotype)
  Trans.dat.2<-cbind(Trans.dat.2, l.phenotype)
}
colnames(Trans.dat.2)<-c("RIX","ID","day","assay_date","IgG1","IgG2ac","IgG2b","IgG3","IgM","TotalG")


for(i in 1:6)
{
  phenotype<-Trans.dat.2[,i+4]
  hist(phenotype, main=abs[i])
}

#D10 transformations
#IgG1 = log
#IgG2ac = sqrt
#IgG2b = sqrt
#IgG3 = log
#IgM = log
#TotalG = none
dat.10=cbind(dat[c(1:4,10)],Trans.dat[c(5,8,9)])
dat.10=cbind(dat.10,Trans.dat.2[c(6,7)])
dat.10=dat.10[col.list]

#visualize final data
for(i in 1:6)
{
  phenotype<-dat.10[,i+4]
  hist(phenotype, main=abs[i],col="lavender")
}

write.csv(dat.10,file=paste("~/Dropbox/Heise/ELISA Antibody/qtls/d10/dat.10.auc.csv",sep=""),row.names=F)

#### day 15 ####

dat=subset(full.dat,full.dat$day==15)

#graph to view data distribution
par(mfrow = c(2,3))

for(i in 1:6)
{
  phenotype<-dat[,i+4]
  phenotype[phenotype<0] <- 0.00000000001
  lm.bc=lm(phenotype~RIX,data=dat)
  bc=boxcox(lm.bc)
  print(bc$x[which.max(bc$y)])
}

for(i in 1:6)
{
  phenotype<-dat[,i+4]
  hist(phenotype, main=abs[i])
}

#see what they look like log transformed - more normal?
Trans.dat<-dat[,1:4]

for(i in 1:6)
{
  phenotype<-dat[,i+4]
  l.phenotype<-log(phenotype,base=10)
  Trans.dat<-cbind(Trans.dat, l.phenotype)
}
colnames(Trans.dat)<-c("RIX","ID","day","assay_date","IgG1","IgG2ac","IgG2b","IgG3","IgM","TotalG")

for(i in 1:6)
{
  phenotype<-Trans.dat[,i+4]
  hist(phenotype, main=abs[i])
}


#see what they look like sqrt transformed - more normal?
Trans.dat.2<-dat[,1:4]

for(i in 1:6)
{
  phenotype<-dat[,i+4]
  l.phenotype<-sqrt(phenotype)
  Trans.dat.2<-cbind(Trans.dat.2, l.phenotype)
}
colnames(Trans.dat.2)<-c("RIX","ID","day","assay_date","IgG1","IgG2ac","IgG2b","IgG3","IgM","TotalG")

for(i in 1:6)
{
  phenotype<-Trans.dat.2[,i+4]
  hist(phenotype, main=abs[i])
}


##D15 Transformations
#IgG1 = sqrt
#IgG2ac = none
#IgG2b = sqrt
#IgG3 = sqrt
#IgM = log
#TotalG = none

dat.15=cbind(dat[c(1:4,6,10)],Trans.dat[c(9)])
dat.15=cbind(dat.15,Trans.dat.2[c(5,7,8)])
dat.15=dat.15[col.list]

#visualize final data
for(i in 1:6)
{
  phenotype<-dat.15[,i+4]
  hist(phenotype, main=abs[i],col="lavender")
}

write.csv(dat.15,paste("~/Dropbox/Heise/ELISA Antibody/qtls/d15/dat.15.auc.csv",sep=""),row.names=F)         


#### day 45 ####
dat=subset(full.dat,full.dat$day==45)

#graph to view data distribution
par(mfrow = c(2,3))

for(i in 1:6)
{
  phenotype<-dat[,i+4]
  phenotype[phenotype<0] <- 0.00000000001
  lm.bc=lm(phenotype~RIX,data=dat)
  bc=boxcox(lm.bc)
  print(bc$x[which.max(bc$y)])
}

for(i in 1:6)
{
  phenotype<-dat[,i+4]
  hist(phenotype, main=abs[i])
}

#see what they look like log transformed - more normal?
Trans.dat<-dat[,1:4]

for(i in 1:6)
{
  phenotype<-dat[,i+4]
  l.phenotype<-log(phenotype, base=10)
  Trans.dat<-cbind(Trans.dat, l.phenotype)
}
colnames(Trans.dat)<-c("RIX","ID","day","assay_date","IgG1","IgG2ac","IgG2b","IgG3","IgM","TotalG")


for(i in 1:6)
{
  phenotype<-Trans.dat[,i+4]
  hist(phenotype, main=abs[i])
}


#see what they look like sqrt transformed - more normal?
Trans.dat.2<-dat[,1:4]

for(i in 1:6)
{
  phenotype<-dat[,i+4]
  l.phenotype<-sqrt(phenotype)
  Trans.dat.2<-cbind(Trans.dat.2, l.phenotype)
}
colnames(Trans.dat.2)<-c("RIX","ID","day","assay_date","IgG1","IgG2ac","IgG2b","IgG3","IgM","TotalG")

for(i in 1:6)
{
  phenotype<-Trans.dat.2[,i+4]
  hist(phenotype, main=abs[i])
}


#see what they look like square transformed - more normal?
Trans.dat.3<-dat[,1:4]

for(i in 1:6)
{
  phenotype<-dat[,i+4]
  l.phenotype<-(phenotype)^2
  Trans.dat.3<-cbind(Trans.dat.3, l.phenotype)
}
colnames(Trans.dat.3)<-c("RIX","ID","day","assay_date","IgG1","IgG2ac","IgG2b","IgG3","IgM","TotalG")

for(i in 1:6)
{
  phenotype<-Trans.dat.3[,i+4]
  hist(phenotype, main=abs[i])
}


#D45 Transformations
#IgG1 = sqrt
#IgG2ac = square
#IgG2b = none
#IgG3 = sqrt
#IgM = log
#TotalG = square
dat.45=cbind(dat[c(1:4,7)],Trans.dat[9])
dat.45=cbind(dat.45,Trans.dat.2[c(5,8)])
dat.45=cbind(dat.45,Trans.dat.3[c(6,10)])
dat.45=dat.45[col.list]

#visualize final data
for(i in 1:6)
{
  phenotype<-dat.45[,i+4]
  hist(phenotype, main=abs[i],col="lavender")
}

write.csv(dat.45,file=paste("~/Dropbox/Heise/ELISA Antibody/qtls/d45/dat.45.auc.csv",sep=""),row.names=F)



##### tests for normality ####
dat=subset(full.dat,full.dat$day==45)
dat=dat.45

for(i in 1:6)
{
  phenotype<-dat[,i+4]
  qqnorm(phenotype, main=abs[i])
}


sha.p=NULL
for(i in 1:6)
{
  phenotype<-dat[,i+4]
  s.p=shapiro.test(phenotype)$p.value
  sha.p <- rbind(sha.p, s.p)
}
rownames(sha.p)=col.list[5:10]
sha.p









###### LAST POSITIVE (>0.2) #######

#load in data to map in wide format
load("wide.lastpos.Rdata")
full.dat<-wide.lastpos
full.dat$RIX=as.factor(full.dat$RIX)

#set list of anitbody isotypes in the column order they appear in
abs=c("IgG1" ,  "IgG2ac" , "IgG2b" , "IgG3" ,  "IgM"  ,  "TotalG")

col.list=c("RIX" , "ID" , "day" , "assay_date", "IgG1" ,"IgG2ac","IgG2b" ,    
           "IgG3","IgM","TotalG")


###### visualize data and transform as necessary (individually by timepoint) #####
## plot data distribution and transform as necessary

#### day 7 ####

dat=subset(full.dat,full.dat$day==7)

# graph to view data distribution #
par(mfrow = c(2,3))

for(i in 1:6)
{
  phenotype<-dat[,i+4]
  hist(phenotype, main=abs[i])
}

#some AUC values are every so slightly below 0.. box cox needs all positive values
#change the values to be positive non-zero numbers so box-cox will run

for(i in 1:6)
{
  phenotype<-dat[,i+4]
  phenotype[phenotype<0] <- 0.00000000001
  phenotype[phenotype==0] <- 0.00000000001
  lm.bc=lm(phenotype~RIX,data=dat)
  bc=boxcox(lm.bc)
  print(bc$x[which.max(bc$y)])
}

#see what they look like log transformed - more normal?
Trans.dat<-dat[,1:4]

for(i in 1:6)
{
  phenotype<-dat[,i+4]
  l.phenotype<-log(phenotype, base=10)
  Trans.dat<-cbind(Trans.dat, l.phenotype)
}
colnames(Trans.dat)<-c("RIX","ID","day","assay_date","IgG1","IgG2ac","IgG2b","IgG3","IgM","TotalG")

for(i in 1:6)
{
  phenotype<-Trans.dat[,i+4]
  hist(phenotype, main=abs[i])
}


#see what they look like sqrt transformed - more normal?
Trans.dat.2<-dat[,1:4]

for(i in 1:6)
{
  phenotype<-dat[,i+4]
  l.phenotype<-sqrt(phenotype)
  Trans.dat.2<-cbind(Trans.dat.2, l.phenotype)
}
colnames(Trans.dat.2)<-c("RIX","ID","day","assay_date","IgG1","IgG2ac","IgG2b","IgG3","IgM","TotalG")

for(i in 1:6)
{
  phenotype<-Trans.dat.2[,i+4]
  hist(phenotype, main=abs[i])
}

##D7 Transformations
#log transform but totalG, which gets sqrt
dat.7=cbind(Trans.dat[c(1:9)],Trans.dat.2[10])

dat.7=dat.7[col.list]

#visualize final data
for(i in 1:6)
{
  phenotype<-dat.7[,i+4]
  hist(phenotype, main=abs[i],col="lavender")
}

write.csv(dat.7,file=paste("~/Dropbox/Heise/ELISA Antibody/qtls/d7/dat.7.lastpos.csv",sep=""),row.names=F)








##

#### day 10 ####
dat=subset(full.dat,full.dat$day==10)

#graph to view data distribution
par(mfrow = c(2,3))

for(i in 1:6)
{
  phenotype<-dat[,i+4]
  hist(phenotype, main=abs[i])
}

for(i in 1:6)
{
  phenotype<-dat[,i+4]
  phenotype[phenotype<0] <- 0.00000000001
  lm.bc=lm(phenotype~RIX,data=dat)
  bc=boxcox(lm.bc)
  print(bc$x[which.max(bc$y)])
}


#see what they look like log transformed - more normal?
Trans.dat<-dat[,1:4]

for(i in 1:6)
{
  phenotype<-dat[,i+4]
  l.phenotype<-log(phenotype, base=10)
  Trans.dat<-cbind(Trans.dat, l.phenotype)
}
colnames(Trans.dat)<-c("RIX","ID","day","assay_date","IgG1","IgG2ac","IgG2b","IgG3","IgM","TotalG")


for(i in 1:6)
{
  phenotype<-Trans.dat[,i+4]
  hist(phenotype, main=abs[i])
}


#see what they look like sqrt transformed - more normal?
Trans.dat.2<-dat[,1:4]

for(i in 1:6)
{
  phenotype<-dat[,i+4]
  l.phenotype<-sqrt(phenotype)
  Trans.dat.2<-cbind(Trans.dat.2, l.phenotype)
}
colnames(Trans.dat.2)<-c("RIX","ID","day","assay_date","IgG1","IgG2ac","IgG2b","IgG3","IgM","TotalG")


for(i in 1:6)
{
  phenotype<-Trans.dat.2[,i+4]
  hist(phenotype, main=abs[i])
}

#D10 transformations
#IgG1 = log
#IgG2ac = sqrt
#IgG2b = sqrt
#IgG3 = log
#IgM = log
#TotalG = none
dat.10=cbind(dat[c(1:4,10)],Trans.dat[c(5,8,9)])
dat.10=cbind(dat.10,Trans.dat.2[c(6,7)])
dat.10=dat.10[col.list]

#visualize final data
for(i in 1:6)
{
  phenotype<-dat.10[,i+4]
  hist(phenotype, main=abs[i],col="lavender")
}

write.csv(dat.10,file=paste("~/Dropbox/Heise/ELISA Antibody/qtls/d10/dat.10.csv",sep=""),row.names=F)

#### day 15 ####

dat=subset(full.dat,full.dat$day==15)

#graph to view data distribution
par(mfrow = c(2,3))

for(i in 1:6)
{
  phenotype<-dat[,i+4]
  phenotype[phenotype<0] <- 0.00000000001
  lm.bc=lm(phenotype~RIX,data=dat)
  bc=boxcox(lm.bc)
  print(bc$x[which.max(bc$y)])
}

for(i in 1:6)
{
  phenotype<-dat[,i+4]
  hist(phenotype, main=abs[i])
}

#see what they look like log transformed - more normal?
Trans.dat<-dat[,1:4]

for(i in 1:6)
{
  phenotype<-dat[,i+4]
  l.phenotype<-log(phenotype,base=10)
  Trans.dat<-cbind(Trans.dat, l.phenotype)
}
colnames(Trans.dat)<-c("RIX","ID","day","assay_date","IgG1","IgG2ac","IgG2b","IgG3","IgM","TotalG")

for(i in 1:6)
{
  phenotype<-Trans.dat[,i+4]
  hist(phenotype, main=abs[i])
}


#see what they look like sqrt transformed - more normal?
Trans.dat.2<-dat[,1:4]

for(i in 1:6)
{
  phenotype<-dat[,i+4]
  l.phenotype<-sqrt(phenotype)
  Trans.dat.2<-cbind(Trans.dat.2, l.phenotype)
}
colnames(Trans.dat.2)<-c("RIX","ID","day","assay_date","IgG1","IgG2ac","IgG2b","IgG3","IgM","TotalG")

for(i in 1:6)
{
  phenotype<-Trans.dat.2[,i+4]
  hist(phenotype, main=abs[i])
}


##D15 Transformations
#IgG1 = sqrt
#IgG2ac = none
#IgG2b = sqrt
#IgG3 = sqrt
#IgM = log
#TotalG = none

dat.15=cbind(dat[c(1:4,6,10)],Trans.dat[c(9)])
dat.15=cbind(dat.15,Trans.dat.2[c(5,7,8)])
dat.15=dat.15[col.list]

#visualize final data
for(i in 1:6)
{
  phenotype<-dat.15[,i+4]
  hist(phenotype, main=abs[i],col="lavender")
}

write.csv(dat.15,paste("~/Dropbox/Heise/ELISA Antibody/qtls/d15/dat.15.csv",sep=""),row.names=F)         


#### day 45 ####
dat=subset(full.dat,full.dat$day==45)

#graph to view data distribution
par(mfrow = c(2,3))

for(i in 1:6)
{
  phenotype<-dat[,i+4]
  phenotype[phenotype<0] <- 0.00000000001
  lm.bc=lm(phenotype~RIX,data=dat)
  bc=boxcox(lm.bc)
  print(bc$x[which.max(bc$y)])
}

for(i in 1:6)
{
  phenotype<-dat[,i+4]
  hist(phenotype, main=abs[i])
}

#see what they look like log transformed - more normal?
Trans.dat<-dat[,1:4]

for(i in 1:6)
{
  phenotype<-dat[,i+4]
  l.phenotype<-log(phenotype, base=10)
  Trans.dat<-cbind(Trans.dat, l.phenotype)
}
colnames(Trans.dat)<-c("RIX","ID","day","assay_date","IgG1","IgG2ac","IgG2b","IgG3","IgM","TotalG")


for(i in 1:6)
{
  phenotype<-Trans.dat[,i+4]
  hist(phenotype, main=abs[i])
}


#see what they look like sqrt transformed - more normal?
Trans.dat.2<-dat[,1:4]

for(i in 1:6)
{
  phenotype<-dat[,i+4]
  l.phenotype<-sqrt(phenotype)
  Trans.dat.2<-cbind(Trans.dat.2, l.phenotype)
}
colnames(Trans.dat.2)<-c("RIX","ID","day","assay_date","IgG1","IgG2ac","IgG2b","IgG3","IgM","TotalG")

for(i in 1:6)
{
  phenotype<-Trans.dat.2[,i+4]
  hist(phenotype, main=abs[i])
}


#see what they look like square transformed - more normal?
Trans.dat.3<-dat[,1:4]

for(i in 1:6)
{
  phenotype<-dat[,i+4]
  l.phenotype<-(phenotype)^2
  Trans.dat.3<-cbind(Trans.dat.3, l.phenotype)
}
colnames(Trans.dat.3)<-c("RIX","ID","day","assay_date","IgG1","IgG2ac","IgG2b","IgG3","IgM","TotalG")

for(i in 1:6)
{
  phenotype<-Trans.dat.3[,i+4]
  hist(phenotype, main=abs[i])
}


#D45 Transformations
#IgG1 = sqrt
#IgG2ac = square
#IgG2b = none
#IgG3 = sqrt
#IgM = log
#TotalG = square
dat.45=cbind(dat[c(1:4,7)],Trans.dat[9])
dat.45=cbind(dat.45,Trans.dat.2[c(5,8)])
dat.45=cbind(dat.45,Trans.dat.3[c(6,10)])
dat.45=dat.45[col.list]

#visualize final data
for(i in 1:6)
{
  phenotype<-dat.45[,i+4]
  hist(phenotype, main=abs[i],col="lavender")
}

write.csv(dat.45,file=paste("~/Dropbox/Heise/ELISA Antibody/qtls/d45/dat.45.auc.csv",sep=""),row.names=F)