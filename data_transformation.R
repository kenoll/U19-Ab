#### setup: load libraries, data format, etc ####
setwd("~/Dropbox/Heise/U19-Ab/")

#load libraries 
library(ggplot2)
library(MASS)

##### AUC #####

#load in data to map in wide format
load("data_bin/wideauc.Rdata")
wideauc$RIX=as.factor(wideauc$RIX)
full.dat<-wideauc

#set list of anitbody isotypes in the column order they appear in
abs=c("IgG1" ,  "IgG2ac" , "IgG2b" , "IgG3" ,  "IgM"  ,  "TotalG")

col.list=colnames(wideauc)

#number of header columns before data
hc=7

###### visualize data and transform as necessary (individually by timepoint) #####
## plot data distribution and transform as necessary

#### day 7 ####

dat=subset(full.dat,full.dat$day==7)

# graph to view data distribution #
par(mfrow = c(2,3))

for(i in 1:6)
{
  phenotype<-dat[,i+hc]
  hist(phenotype, main=abs[i])
}

#some AUC values are every so slightly below 0.. box cox needs all positive values
#change the values to be positive non-zero numbers so box-cox will run

for(i in 1:6)
{
  phenotype<-dat[,i+hc]
  phenotype[phenotype<0] <- 0.00000000001
  lm.bc=lm(phenotype~RIX,data=dat)
  bc=boxcox(lm.bc)
  print(bc$x[which.max(bc$y)])
}

#see what they look like log transformed - more normal?
dat.log<-dat[,1:hc]

for(i in 1:6)
{
  phenotype<-dat[,i+hc]
  l.phenotype<-log(phenotype, base=10)
  dat.log<-cbind(dat.log, l.phenotype)
}
colnames(dat.log)<-col.list

for(i in 1:6)
{
  phenotype<-dat.log[,i+hc]
  hist(phenotype, main=abs[i])
}


#see what they look like sqrt transformed - more normal?
dat.sqrt<-dat[,1:hc]

for(i in 1:6)
{
  phenotype<-dat[,i+hc]
  l.phenotype<-sqrt(phenotype)
  dat.sqrt<-cbind(dat.sqrt, l.phenotype)
}
colnames(dat.sqrt)<-col.list

for(i in 1:6)
{
  phenotype<-dat.sqrt[,i+hc]
  hist(phenotype, main=abs[i])
}


##D7 Transformations
#log transform all
dat.7=cbind(dat.log[c(1:13)])
dat.7=dat.7[col.list]

#visualize final data
for(i in 1:6)
{
  phenotype<-dat.7[,i+hc]
  hist(phenotype, main=abs[i],col="lavender")
}

write.csv(dat.7,file=paste("~/Dropbox/Heise/U19-Ab/data_bin/d7/dat.7.auc.csv",sep=""),row.names=F)


#### day 10 ####
dat=subset(full.dat,full.dat$day==10)

#graph to view data distribution
par(mfrow = c(2,3))

for(i in 1:6)
{
  phenotype<-dat[,i+hc]
  hist(phenotype, main=abs[i])
}

for(i in 1:6)
{
  phenotype<-dat[,i+hc]
  phenotype[phenotype<0] <- 0.00000000001
  lm.bc=lm(phenotype~RIX,data=dat)
  bc=boxcox(lm.bc)
  print(bc$x[which.max(bc$y)])
}


#see what they look like log transformed - more normal?
dat.log<-dat[,1:hc]

for(i in 1:6)
{
  phenotype<-dat[,i+hc]
  l.phenotype<-log(phenotype, base=10)
  dat.log<-cbind(dat.log, l.phenotype)
}
colnames(dat.log)<-col.list


for(i in 1:6)
{
  phenotype<-dat.log[,i+hc]
  hist(phenotype, main=abs[i])
}


#see what they look like sqrt transformed - more normal?
dat.sqrt<-dat[,1:hc]

for(i in 1:6)
{
  phenotype<-dat[,i+hc]
  l.phenotype<-sqrt(phenotype)
  dat.sqrt<-cbind(dat.sqrt, l.phenotype)
}
colnames(dat.sqrt)<-col.list


for(i in 1:6)
{
  phenotype<-dat.sqrt[,i+hc]
  hist(phenotype, main=abs[i])
}

#D10 transformations
#IgG1 = log
#IgG2ac = sqrt
#IgG2b = sqrt
#IgG3 = log
#IgM = log
#TotalG = none
dat.10=dat[1:hc]
dat.10=cbind(dat.10,dat["TotalG"])  
dat.10=cbind(dat.10,dat.log[c("IgG1","IgG3","IgM")])
dat.10=cbind(dat.10,dat.sqrt[c("IgG2ac","IgG2b")])
dat.10=dat.10[col.list]

#visualize final data
for(i in 1:6)
{
  phenotype<-dat.10[,i+hc]
  hist(phenotype, main=abs[i],col="lavender")
}

write.csv(dat.10,file=paste("~/Dropbox/Heise/U19-Ab/data_bin/d10/dat.10.auc.csv",sep=""),row.names=F)

#### day 15 ####

dat=subset(full.dat,full.dat$day==15)

#graph to view data distribution
par(mfrow = c(2,3))

for(i in 1:6)
{
  phenotype<-dat[,i+hc]
  phenotype[phenotype<0] <- 0.00000000001
  lm.bc=lm(phenotype~RIX,data=dat)
  bc=boxcox(lm.bc)
  print(bc$x[which.max(bc$y)])
}

for(i in 1:6)
{
  phenotype<-dat[,i+hc]
  hist(phenotype, main=abs[i])
}

#see what they look like log transformed - more normal?
dat.log<-dat[,1:hc]

for(i in 1:6)
{
  phenotype<-dat[,i+hc]
  l.phenotype<-log(phenotype,base=10)
  dat.log<-cbind(dat.log, l.phenotype)
}
colnames(dat.log)=col.list

for(i in 1:6)
{
  phenotype<-dat.log[,i+hc]
  hist(phenotype, main=abs[i])
}


#see what they look like sqrt transformed - more normal?
dat.sqrt<-dat[,1:hc]

for(i in 1:6)
{
  phenotype<-dat[,i+hc]
  l.phenotype<-sqrt(phenotype)
  dat.sqrt<-cbind(dat.sqrt, l.phenotype)
}
colnames(dat.sqrt)=col.list

for(i in 1:6)
{
  phenotype<-dat.sqrt[,i+hc]
  hist(phenotype, main=abs[i])
}


##D15 Transformations
#IgG1 = sqrt
#IgG2ac = none
#IgG2b = sqrt
#IgG3 = sqrt
#IgM = log
#TotalG = none

dat.15=dat[1:hc]
dat.15=cbind(dat.15,dat[c("IgG2ac","TotalG")])
dat.15=cbind(dat.15,dat.log["IgM"])
dat.15=cbind(dat.15,dat.sqrt[c("IgG1","IgG2b","IgG3")])
dat.15=dat.15[col.list]

#visualize final data
for(i in 1:6)
{
  phenotype<-dat.15[,i+hc]
  hist(phenotype, main=abs[i],col="lavender")
}

write.csv(dat.15,paste("~/Dropbox/Heise/U19-Ab/data_bin/d15/dat.15.auc.csv",sep=""),row.names=F)         


#### day 45 ####
dat=subset(full.dat,full.dat$day==45)

#graph to view data distribution
par(mfrow = c(2,3))

for(i in 1:6)
{
  phenotype<-dat[,i+hc]
  phenotype[phenotype<0] <- 0.00000000001
  lm.bc=lm(phenotype~RIX,data=dat)
  bc=boxcox(lm.bc)
  print(bc$x[which.max(bc$y)])
}

for(i in 1:6)
{
  phenotype<-dat[,i+hc]
  hist(phenotype, main=abs[i])
}

#see what they look like log transformed - more normal?
dat.log<-dat[,1:hc]

for(i in 1:6)
{
  phenotype<-dat[,i+hc]
  l.phenotype<-log(phenotype, base=10)
  dat.log<-cbind(dat.log, l.phenotype)
}
colnames(dat.log)=col.list


for(i in 1:6)
{
  phenotype<-dat.log[,i+hc]
  hist(phenotype, main=abs[i])
}


#see what they look like sqrt transformed - more normal?
dat.sqrt<-dat[,1:hc]

for(i in 1:6)
{
  phenotype<-dat[,i+hc]
  l.phenotype<-sqrt(phenotype)
  dat.sqrt<-cbind(dat.sqrt, l.phenotype)
}
colnames(dat.sqrt)=col.list

for(i in 1:6)
{
  phenotype<-dat.sqrt[,i+hc]
  hist(phenotype, main=abs[i])
}


#see what they look like square transformed - more normal?
dat.sq2<-dat[,1:hc]

for(i in 1:6)
{
  phenotype<-dat[,i+hc]
  l.phenotype<-(phenotype)^2
  dat.sq2<-cbind(dat.sq2, l.phenotype)
}
colnames(dat.sq2)=col.list

for(i in 1:6)
{
  phenotype<-dat.sq2[,i+hc]
  hist(phenotype, main=abs[i])
}


#D45 Transformations
#IgG1 = sqrt
#IgG2ac = square
#IgG2b = none
#IgG3 = sqrt
#IgM = log
#TotalG = square
dat.45=dat[1:hc]
dat.45=cbind(dat.45,dat["IgG2b"])
dat.45=cbind(dat.45,dat.log["IgM"])
dat.45=cbind(dat.45,dat.sqrt[c("IgG1","IgG3")])
dat.45=cbind(dat.45,dat.sq2[c("IgG2ac","TotalG")])
dat.45=dat.45[col.list]

#visualize final data
for(i in 1:6)
{
  phenotype<-dat.45[,i+hc]
  hist(phenotype, main=abs[i],col="lavender")
}

write.csv(dat.45,file=paste("~/Dropbox/Heise/U19-Ab/data_bin/d45/dat.45.auc.csv",sep=""),row.names=F)



##### tests for normality ####
dat=subset(full.dat,full.dat$day==45)
dat=dat.45

for(i in 1:6)
{
  phenotype<-dat[,i+hc]
  qqnorm(phenotype, main=abs[i])
}


sha.p=NULL
for(i in 1:6)
{
  phenotype<-dat[,i+hc]
  s.p=shapiro.test(phenotype)$p.value
  sha.p <- rbind(sha.p, s.p)
}
rownames(sha.p)=col.list[(hc+1):ncol(dat)]
sha.p



###### LAST POSITIVE (>0.2) #######

#load in data to map in wide format
load("data_bin/wide.lastpos.Rdata")
full.dat<-wide.lastpos
full.dat$RIX=as.factor(full.dat$RIX)

#set list of anitbody isotypes in the column order they appear in
abs=c("IgG1" ,  "IgG2ac" , "IgG2b" , "IgG3" ,  "IgM"  ,  "TotalG")

col.list=colnames(full.dat)


###### visualize data and transform as necessary (individually by timepoint) #####
## plot data distribution and transform as necessary

#### day 7 ####

dat=subset(full.dat,full.dat$day==7)

# graph to view data distribution #
par(mfrow = c(2,3))

for(i in 1:6)
{
  phenotype<-dat[,i+hc]
  hist(phenotype, main=abs[i])
}

#some AUC values are every so slightly below 0.. box cox needs all positive values
#change the values to be positive non-zero numbers so box-cox will run

for(i in 1:6)
{
  phenotype<-dat[,i+hc]
  phenotype[phenotype<0] <- 0.00000000001
  phenotype[phenotype==0] <- 0.00000000001
  lm.bc=lm(phenotype~RIX,data=dat)
  bc=boxcox(lm.bc)
  print(bc$x[which.max(bc$y)])
}

#see what they look like log transformed - more normal?
dat.log<-dat[,1:hc]

for(i in 1:6)
{
  phenotype<-dat[,i+hc]
  l.phenotype<-log(phenotype, base=10)
  dat.log<-cbind(dat.log, l.phenotype)
}
colnames(dat.log)=col.list

for(i in 1:6)
{
  phenotype<-dat.log[,i+hc]
  hist(phenotype, main=abs[i])
}


#see what they look like sqrt transformed - more normal?
dat.sqrt<-dat[,1:hc]

for(i in 1:6)
{
  phenotype<-dat[,i+hc]
  l.phenotype<-sqrt(phenotype)
  dat.sqrt<-cbind(dat.sqrt, l.phenotype)
}
colnames(dat.sqrt)=col.list

for(i in 1:6)
{
  phenotype<-dat.sqrt[,i+hc]
  hist(phenotype, main=abs[i])
}

##D7 Transformations
#log transform but totalG, which gets sqrt
dat.7=dat[1:hc]
dat.7=cbind(dat.7,dat.log[c("IgG1","IgG2ac","IgG2b","IgG3","IgM")])
dat.7=cbind(dat.7,dat.sqrt["TotalG"])
dat.7=dat.7[col.list]

#visualize final data
for(i in 1:6)
{
  phenotype<-dat.7[,i+hc]
  hist(phenotype, main=abs[i],col="lavender")
}

write.csv(dat.7,file=paste("~/Dropbox/Heise/U19-Ab/data_bin/d7/dat.7.lastpos.csv",sep=""),row.names=F)


##

#### day 10 ####
dat=subset(full.dat,full.dat$day==10)

#graph to view data distribution
par(mfrow = c(2,3))

for(i in 1:6)
{
  phenotype<-dat[,i+hc]
  hist(phenotype, main=abs[i])
}

for(i in 1:6)
{
  phenotype<-dat[,i+hc]
  phenotype[phenotype<0] <- 0.00000000001
  phenotype[phenotype==0] <- 0.00000000001
  lm.bc=lm(phenotype~RIX,data=dat)
  bc=boxcox(lm.bc)
  print(bc$x[which.max(bc$y)])
}


#see what they look like log transformed - more normal?
dat.log<-dat[,1:hc]

for(i in 1:6)
{
  phenotype<-dat[,i+hc]
  l.phenotype<-log(phenotype, base=10)
  dat.log<-cbind(dat.log, l.phenotype)
}
colnames(dat.log)=col.list


for(i in 1:6)
{
  phenotype<-dat.log[,i+hc]
  hist(phenotype, main=abs[i])
}


#see what they look like sqrt transformed - more normal?
dat.sqrt<-dat[,1:hc]

for(i in 1:6)
{
  phenotype<-dat[,i+hc]
  l.phenotype<-sqrt(phenotype)
  dat.sqrt<-cbind(dat.sqrt, l.phenotype)
}
colnames(dat.sqrt)=col.list


for(i in 1:6)
{
  phenotype<-dat.sqrt[,i+hc]
  hist(phenotype, main=abs[i])
}

#D10 transformations
#IgG1 = log
#IgG2ac = sqrt
#IgG2b = sqrt
#IgG3 = log
#IgM = log
#TotalG = none
dat.10=dat[1:hc]
dat.10=cbind(dat.10,dat["TotalG"])  
dat.10=cbind(dat.10,dat.log[c("IgG1","IgG3","IgM")])
dat.10=cbind(dat.10,dat.sqrt[c("IgG2ac","IgG2b")])
dat.10=dat.10[col.list]

#visualize final data
for(i in 1:6)
{
  phenotype<-dat.10[,i+hc]
  hist(phenotype, main=abs[i],col="lavender")
}

write.csv(dat.10,file=paste("~/Dropbox/Heise/U19-Ab/data_bin/d10/dat.10.lastpos.csv",sep=""),row.names=F)

#### day 15 ####

dat=subset(full.dat,full.dat$day==15)

#graph to view data distribution
par(mfrow = c(2,3))

for(i in 1:6)
{
  phenotype<-dat[,i+hc]
  phenotype[phenotype<0] <- 0.00000000001
  phenotype[phenotype==0] <- 0.00000000001
  lm.bc=lm(phenotype~RIX,data=dat)
  bc=boxcox(lm.bc)
  print(bc$x[which.max(bc$y)])
}

for(i in 1:6)
{
  phenotype<-dat[,i+hc]
  hist(phenotype, main=abs[i])
}

#see what they look like log transformed - more normal?
dat.log<-dat[,1:hc]

for(i in 1:6)
{
  phenotype<-dat[,i+hc]
  l.phenotype<-log(phenotype,base=10)
  dat.log<-cbind(dat.log, l.phenotype)
}
colnames(dat.log)=col.list

for(i in 1:6)
{
  phenotype<-dat.log[,i+hc]
  hist(phenotype, main=abs[i])
}


#see what they look like sqrt transformed - more normal?
dat.sqrt<-dat[,1:hc]

for(i in 1:6)
{
  phenotype<-dat[,i+hc]
  l.phenotype<-sqrt(phenotype)
  dat.sqrt<-cbind(dat.sqrt, l.phenotype)
}
colnames(dat.sqrt)=col.list

for(i in 1:6)
{
  phenotype<-dat.sqrt[,i+hc]
  hist(phenotype, main=abs[i])
}


#see what they look like sq^2 transformed - more normal?
dat.sq2<-dat[,1:hc]

for(i in 1:6)
{
  phenotype<-dat[,i+hc]
  l.phenotype<-(phenotype)^2
  dat.sq2<-cbind(dat.sq2, l.phenotype)
}
colnames(dat.sq2)=col.list

for(i in 1:6)
{
  phenotype<-dat.sqrt[,i+hc]
  hist(phenotype, main=abs[i])
}




##D15 Transformations
#IgG1 = sqrt
#IgG2ac = square
#IgG2b = none
#IgG3 = sqrt
#IgM = log
#TotalG = square

dat.15=dat[1:hc]
dat.15=cbind(dat.15,dat[c("IgG2b")])
dat.15=cbind(dat.15,dat.log["IgM"])
dat.15=cbind(dat.15,dat.sqrt[c("IgG1","IgG3")])
dat.15=cbind(dat.15,dat.sq2[c("IgG2ac","TotalG")])
dat.15=dat.15[col.list]

#visualize final data
for(i in 1:6)
{
  phenotype<-dat.15[,i+hc]
  hist(phenotype, main=abs[i],col="lavender")
}

write.csv(dat.15,paste("~/Dropbox/Heise/U19-Ab/data_bin/d15/dat.15.lastpos.csv",sep=""),row.names=F)         


#### day 45 ####
dat=subset(full.dat,full.dat$day==45)

#graph to view data distribution
par(mfrow = c(2,3))

for(i in 1:6)
{
  phenotype<-dat[,i+hc]
  phenotype[phenotype<0] <- 0.00000000001
  phenotype[phenotype==0] <- 0.00000000001
  lm.bc=lm(phenotype~RIX,data=dat)
  bc=boxcox(lm.bc)
  print(bc$x[which.max(bc$y)])
}

for(i in 1:6)
{
  phenotype<-dat[,i+hc]
  hist(phenotype, main=abs[i])
}

#see what they look like log transformed - more normal?
dat.log<-dat[,1:hc]

for(i in 1:6)
{
  phenotype<-dat[,i+hc]
  l.phenotype<-log(phenotype, base=10)
  dat.log<-cbind(dat.log, l.phenotype)
}
colnames(dat.log)=col.list


for(i in 1:6)
{
  phenotype<-dat.log[,i+hc]
  hist(phenotype, main=abs[i])
}


#see what they look like sqrt transformed - more normal?
dat.sqrt<-dat[,1:hc]

for(i in 1:6)
{
  phenotype<-dat[,i+hc]
  l.phenotype<-sqrt(phenotype)
  dat.sqrt<-cbind(dat.sqrt, l.phenotype)
}
colnames(dat.sqrt)=col.list

for(i in 1:6)
{
  phenotype<-dat.sqrt[,i+hc]
  hist(phenotype, main=abs[i])
}


#see what they look like square transformed - more normal?
dat.sq2<-dat[,1:hc]

for(i in 1:6)
{
  phenotype<-dat[,i+hc]
  l.phenotype<-(phenotype)^2
  dat.sq2<-cbind(dat.sq2, l.phenotype)
}
colnames(dat.sq2)=col.list

for(i in 1:6)
{
  phenotype<-dat.sq2[,i+hc]
  hist(phenotype, main=abs[i])
}


#D45 Transformations
#IgG1 = sqrt
#IgG2ac = square
#IgG2b = square
#IgG3 = sqrt
#IgM = log
#TotalG = square
dat.45=dat[1:hc]
dat.45=cbind(dat.45,dat.log["IgM"])
dat.45=cbind(dat.45,dat.sqrt[c("IgG1","IgG3")])
dat.45=cbind(dat.45,dat.sq2[c("IgG2ac","IgG2b","TotalG")])
dat.45=dat.45[col.list]

#visualize final data
for(i in 1:6)
{
  phenotype<-dat.45[,i+hc]
  hist(phenotype, main=abs[i],col="lavender")
}

write.csv(dat.45,file=paste("~/Dropbox/Heise/U19-Ab/data_bin/d45/dat.45.lastpos.csv",sep=""),row.names=F)
