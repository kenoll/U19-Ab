#### setup: load libraries, data format, etc ####
setwd("~/Dropbox/Heise/U19-Ab/sars_ab/")

#load libraries 
library(ggplot2)
library(MASS)

##### AUC #####

#load in data to map in wide format
wideauc=read.csv("data_bin/AUC_by_isotype.csv")
wideauc$RIX=as.factor(wideauc$RIX)
full.dat<-wideauc
full.dat$day = full.dat$day %>% relevel("DRECHAL") %>% relevel("D29") %>% relevel("D15") %>% relevel("D10") %>% relevel("D7")

#set list of anitbody isotypes in the column order they appear in
abs=c("IgG1" ,  "IgG2ac" , "IgG2b" , "IgG3" ,  "IgM"  ,  "TotalG")

col.list=colnames(wideauc)

#number of header columns before data
hc=6

###### visualize data and transform as necessary (individually by timepoint) #####
## plot data distribution and transform as necessary

#### day 7 ####

dat=subset(full.dat,full.dat$day=="D7")
apply(dat[,6:11],2,range,na.rm=T) #see if you can get away with making all very close to zero numbers zero

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
  phenotype[phenotype<0] <- 0.00000000001
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
  phenotype[phenotype<0] <- 0.00000000001
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
dat.7=cbind(dat.log[c(1:12)])
dat.7=dat.7[col.list]

#visualize final data
for(i in 1:6)
{
  phenotype<-dat.7[,i+hc]
  hist(phenotype, main=abs[i],col="lavender")
}

write.csv(dat.7,file=paste("data_bin/dat.7.auc.100.csv",sep=""),row.names=F)


#### day 10 ####
dat=subset(full.dat,full.dat$day=="D10")

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
  phenotype[phenotype<0] <- 0.00000000001
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
  phenotype[phenotype<0] <- 0.00000000001
  l.phenotype<-sqrt(phenotype)
  dat.sqrt<-cbind(dat.sqrt, l.phenotype)
}
colnames(dat.sqrt)<-col.list


for(i in 1:6)
{
  phenotype<-dat.sqrt[,i+hc]
  hist(phenotype, main=abs[i])
}


#see what they look like inversed - more normal?
dat.inv<-dat[,1:hc]

for(i in 1:6)
{
  phenotype<-dat[,i+hc]
  phenotype[phenotype<0] <- 0.00000000001
  l.phenotype<-(phenotype)^-1
  dat.inv<-cbind(dat.inv, l.phenotype)
}
colnames(dat.inv)<-col.list


for(i in 1:6)
{
  phenotype<-dat.inv[,i+hc]
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

write.csv(dat.10,file=paste("data_bin/dat.10.auc.100.csv",sep=""),row.names=F)

#### day 15 ####

dat=subset(full.dat,full.dat$day=="D15")

#graph to view data distribution
par(mfrow = c(2,2))

for(i in c(1,3,4,5))
{
  phenotype<-dat[,i+hc]
  phenotype[phenotype<0] <- 0.00000000001
  lm.bc=lm(phenotype~RIX,data=dat)
  bc=boxcox(lm.bc)
  print(bc$x[which.max(bc$y)])
}

for(i in c(1,3,4,5))
{
  phenotype<-dat[,i+hc]
  hist(phenotype, main=abs[i])
}

#see what they look like log transformed - more normal?
dat.log<-dat[,1:hc]

for(i in c(1,3,4,5))
{
  phenotype<-dat[,i+hc]
  phenotype[phenotype<0] <- 0.00000000001
  l.phenotype<-log(phenotype,base=10)
  dat.log<-cbind(dat.log, l.phenotype)
}
colnames(dat.log)=col.list[-c(hc+2,hc+6)]

for(i in c(1:4))
{
  phenotype<-dat.log[,i+hc]
  hist(phenotype, main=abs[i])
}


#see what they look like sqrt transformed - more normal?
dat.sqrt<-dat[,1:hc]

for(i in c(1,3,4,5))
{
  phenotype<-dat[,i+hc]
  phenotype[phenotype<0] <- 0.00000000001
  l.phenotype<-sqrt(phenotype)
  dat.sqrt<-cbind(dat.sqrt, l.phenotype)
}
colnames(dat.sqrt)=col.list[-c(hc+2,hc+6)]

for(i in 1:4)
{
  phenotype<-dat.sqrt[,i+hc]
  hist(phenotype, main=abs[i])
}


##D15 Transformations
#IgG1 = sqrt
#IgG2ac = NA
#IgG2b = sqrt
#IgG3 = sqrt
#IgM = log
#TotalG = NA

dat.15=dat[1:hc]
dat.15=cbind(dat.15,dat.log["IgM"])
dat.15=cbind(dat.15,dat.sqrt[c("IgG1","IgG2b","IgG3")])
dat.15=dat.15[col.list[-c(hc+2,hc+6)]]

#visualize final data
for(i in 1:4)
{
  phenotype<-dat.15[,i+hc]
  hist(phenotype, main=abs[i],col="lavender")
}

write.csv(dat.15,paste("data_bin/dat.15.auc.100.csv",sep=""),row.names=F)         


#### day 29 ####
dat=subset(full.dat,full.dat$day=="D29")

#graph to view data distribution
par(mfrow = c(1))

for(i in 1:1)
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

for(i in 5)
{
  phenotype<-dat[,i+hc]
  phenotype[phenotype<0] <- 0.00000000001
  l.phenotype<-log(phenotype, base=10)
  dat.log<-cbind(dat.log, l.phenotype)
}
colnames(dat.log)=col.list[c(1:hc,hc+5)]


for(i in 1:6)
{
  phenotype<-dat.log[,i+hc]
  hist(phenotype, main=abs[i])
}


#see what they look like sqrt transformed - more normal?
dat.sqrt<-dat[,1:hc]

for(i in 5)
{
  phenotype<-dat[,i+hc]
  phenotype[phenotype<0] <- 0.00000000001
  l.phenotype<-sqrt(phenotype)
  dat.sqrt<-cbind(dat.sqrt, l.phenotype)
}
colnames(dat.sqrt)=col.list[c(1:hc,hc+5)]

for(i in 1:1)
{
  phenotype<-dat.sqrt[,i+hc]
  hist(phenotype, main="abs[i]")
}


# #see what they look like square transformed - more normal?
# dat.sq2<-dat[,1:hc]
# 
# for(i in 1:6)
# {
#   phenotype<-dat[,i+hc]
#   l.phenotype<-(phenotype)^2
#   dat.sq2<-cbind(dat.sq2, l.phenotype)
# }
# colnames(dat.sq2)=col.list
# 
# for(i in 1:6)
# {
#   phenotype<-dat.sq2[,i+hc]
#   hist(phenotype, main=abs[i])
# }


#D45 Transformations
#IgG1 = NA
#IgG2ac = NA
#IgG2b = NA
#IgG3 = NA
#IgM = log
#TotalG = NA
dat.29=dat[1:hc]
dat.29=cbind(dat.29,dat.log["IgM"])
dat.29=dat.29[col.list[c(1:hc,hc+5)]]


#visualize final data
for(i in 1:1)
{
  phenotype<-dat.29[,i+hc]
  hist(phenotype, main=abs[i],col="lavender")
}

write.csv(dat.29,file=paste("data_bin/dat.29.auc.100.csv",sep=""),row.names=F)


#### rechallenge ####
dat=subset(full.dat,full.dat$day=="DRECHAL")

#graph to view data distribution
par(mfrow = c(1,1))

for(i in 1:1)
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

for(i in 5)
{
  phenotype<-dat[,i+hc]
  phenotype[phenotype<0] <- 0.00000000001
  l.phenotype<-log(phenotype, base=10)
  dat.log<-cbind(dat.log, l.phenotype)
}
colnames(dat.log)=col.list[c(1:hc,hc+5)]


for(i in 1:6)
{
  phenotype<-dat.log[,i+hc]
  hist(phenotype, main=abs[i])
}


#see what they look like sqrt transformed - more normal?
dat.sqrt<-dat[,1:hc]

for(i in 5)
{
  phenotype<-dat[,i+hc]
  phenotype[phenotype<0] <- 0.00000000001
  l.phenotype<-sqrt(phenotype)
  dat.sqrt<-cbind(dat.sqrt, l.phenotype)
}
colnames(dat.sqrt)=col.list[c(1:hc,hc+5)]

for(i in 1:1)
{
  phenotype<-dat.sqrt[,i+hc]
  hist(phenotype, main="abs[i]")
}


# #see what they look like square transformed - more normal?
# dat.sq2<-dat[,1:hc]
# 
# for(i in 1:6)
# {
#   phenotype<-dat[,i+hc]
#   l.phenotype<-(phenotype)^2
#   dat.sq2<-cbind(dat.sq2, l.phenotype)
# }
# colnames(dat.sq2)=col.list
# 
# for(i in 1:6)
# {
#   phenotype<-dat.sq2[,i+hc]
#   hist(phenotype, main=abs[i])
# }


#D45 Transformations
#IgG1 = NA
#IgG2ac = NA
#IgG2b = NA
#IgG3 = NA
#IgM = log
#TotalG = NA
dat.rechal=dat[1:hc]
dat.rechal=cbind(dat.rechal,dat.log["IgM"])
dat.rechal=dat.rechal[col.list[c(1:hc,hc+5)]]



##### tests for normality ####
dat=subset(full.dat,full.dat$day=="D7")
dat=dat.7

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


