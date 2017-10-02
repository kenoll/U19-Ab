#### setup: load libraries, data format, etc ####
setwd("~/Dropbox/Heise/U19-Ab/antibody")

#load libraries 
library(ggplot2)
library(doBy)
library(reshape2)
library(MASS)
library(dplyr)

#load in wide format data to map
dat=read.csv("data_bin/AUC_by_isotype.csv")

#set list of anitbody isotypes in the column order they appear in
abs=c("IgG1" ,"IgG2ac", "IgG2b", "IgG3", "IgM", "TotalG")
col.list=colnames(dat)

#number of header columns before data
hc=10

kin=summaryBy(IgG1+IgG2ac+IgG2b+IgG3+IgM+TotalG ~ 
                     RIX + day, data=dat, FUN=mean, na.rm=T)
colnames(kin)[3:8] = gsub(".mean","",colnames(kin[3:8]))

kin=melt(kin,id.vars=c("RIX","day"),variable.name="isotype",value.name="AUC")
kin=dcast(kin, RIX + isotype ~ day,mean,value.var='AUC')
colnames(kin)[3:6]=c("d7","d10","d15","d45")

kin.method="subtract"

kin$d7_d10=kin$d10-kin$d7
kin$d10_d15=kin$d15-kin$d10
kin$d15_d45=kin$d45-kin$d15

kin=kin[c(1:2,7:9)]

kin=melt(kin,id.vars=c("RIX","isotype"),measure.vars=c("d7_d10","d10_d15","d15_d45"),variable.name="day",value.name="AUC")
kin=dcast(kin, RIX + day ~ isotype,mean,value.var='AUC')

#add minimum value so all values positive
kin[3:8]=kin[3:8]+(abs(min(kin[3:8],na.rm=T))+0.01)

write.csv(kin,paste0("phenotypes/kinetics/auc_kinetics_",kin.method,"_raw.csv"))


#######

full.dat=kin

#set list of anitbody isotypes in the column order they appear in
abs=c("IgG1" ,"IgG2ac", "IgG2b", "IgG3", "IgM", "TotalG")
col.list=colnames(full.dat)

#number of header columns before data
hc=2

#### day 7 to 10  ####
dat=subset(full.dat,full.dat$day=="d7_d10")

# baseline distrubtion
par(mfrow = c(2,3))
hist.ab(dat,hc)

#some AUC values are every so slightly below 0.. box cox needs all positive values
#change the values to be positive non-zero numbers so box-cox will run
# box.ab(dat,hc=hc,num=6)
# for(i in 1:6)
# {
#   phenotype<-dat[,i+hc]
#   lm.bc=lm(phenotype~RIX,data=dat)
#   bc=boxcox(lm.bc)
#   print(bc$x[which.max(bc$y)])
# }

#test various transformations
dat.inv2=trans.ab(dat,-2,hc)
  colnames(dat.inv2)<-col.list
  hist.ab(dat.inv2,hc)  

dat.inv=trans.ab(dat,-1,hc)
  colnames(dat.inv)<-col.list
  hist.ab(dat.inv,hc)  

dat.invsqrt=trans.ab(dat,-0.5,hc)
  colnames(dat.invsqrt)<-col.list
  hist.ab(dat.invsqrt,hc)  
  
dat.log<-log.ab(dat,hc)
  colnames(dat.log)<-col.list
  hist.ab(dat.log,hc)

dat.sqrt=trans.ab(dat,0.5,hc)
  colnames(dat.sqrt)<-col.list
  hist.ab(dat.sqrt,hc)

dat.sq2=trans.ab(dat,2,hc)
  colnames(dat.sq2)<-col.list
  hist.ab(dat.sq2,hc)
  

##d7-10 Transformations
# #for ratio - log transform all
# dat.7_10=dat[1:hc] %>%
#   cbind(dat.log[c("IgG1","IgG2ac","IgG2b","IgG3","IgM","TotalG")])
# dat.7_10=dat.7_10[col.list]

#for subtration
dat.7_10=dat[1:hc] %>%
  cbind(dat.log[c("IgG1","IgG3","IgM","TotalG")]) %>%
  cbind(dat.inv["IgG2ac"]) %>%
  cbind(dat.inv2["IgG2b"])
dat.7_10=dat.7_10[col.list]

#visualize final data
par(mfrow=c(2,3))
for(i in 1:6)
{phenotype<-dat.7_10[,i+hc];hist(phenotype, main=paste0("Kinetics: Day 7 to 10 ",abs[i]),col="honeydew")}
write.csv(dat.7_10,file=paste0("data_bin/kinetics/dat.7-10_",kin.method,"_auc.csv"),row.names=F)


#### day 10 to 15 ####
dat=subset(full.dat,full.dat$day=="d10_d15")

# baseline distrubtion
par(mfrow = c(2,3))
hist.ab(dat,hc)

#test various transformations
dat.inv2=trans.ab(dat,-2,hc)
colnames(dat.inv2)<-col.list
hist.ab(dat.inv2,hc)  

dat.inv=trans.ab(dat,-1,hc)
colnames(dat.inv)<-col.list
hist.ab(dat.inv,hc)  

dat.invsqrt=trans.ab(dat,-0.5,hc)
colnames(dat.invsqrt)<-col.list
hist.ab(dat.invsqrt,hc)  

dat.log<-log.ab(dat,hc)
colnames(dat.log)<-col.list
hist.ab(dat.log,hc)

dat.sqrt=trans.ab(dat,0.5,hc)
colnames(dat.sqrt)<-col.list
hist.ab(dat.sqrt,hc)

dat.sq2=trans.ab(dat,2,hc)
colnames(dat.sq2)<-col.list
hist.ab(dat.sq2,hc)


##d10_d15 Transformations
# #ratio: log transform all
# dat.10_15=dat[1:hc] %>%
#   cbind(dat.log[c("IgG1","IgG2ac","IgG2b","IgG3","IgM","TotalG")])
# dat.10_15=dat.10_15[col.list]

#subtraction
dat.10_15=dat[1:hc] %>%
cbind(dat[c("IgG2ac","IgG2b","IgG3","IgM","TotalG")]) %>%
  cbind(dat.invsqrt["IgG1"])
  dat.10_15=dat.10_15[col.list]

#visualize final data
par(mfrow=c(2,3))
for(i in 1:6)
{phenotype<-dat.10_15[,i+hc];hist(phenotype, main=paste0("Kinetics: Day 10 to 15 ",abs[i]),col="honeydew")}
write.csv(dat.10_15,file=paste0("data_bin/kinetics/dat.10-15_",kin.method,"_auc.csv"),row.names=F)


#### day 15 to 45 ####
dat=subset(full.dat,full.dat$day=="d15_d45")

# baseline distrubtion
par(mfrow = c(2,3))
hist.ab(dat,hc)

#test various transformations
dat.inv2=trans.ab(dat,-2,hc)
colnames(dat.inv2)<-col.list
hist.ab(dat.inv2,hc)  

dat.inv=trans.ab(dat,-1,hc)
colnames(dat.inv)<-col.list
hist.ab(dat.inv,hc)  

dat.invsqrt=trans.ab(dat,-0.5,hc)
colnames(dat.invsqrt)<-col.list
hist.ab(dat.invsqrt,hc)  

dat.log<-log.ab(dat,hc)
colnames(dat.log)<-col.list
hist.ab(dat.log,hc)

dat.sqrt=trans.ab(dat,0.5,hc)
colnames(dat.sqrt)<-col.list
hist.ab(dat.sqrt,hc)

dat.sq2=trans.ab(dat,2,hc)
colnames(dat.sq2)<-col.list
hist.ab(dat.sq2,hc)

##d15_d45 Transformations
# #ratio: log transform all
# dat.15_45=dat[1:hc] %>%
#   cbind(dat.log[c("IgG1","IgG2ac","IgG2b","IgG3","IgM","TotalG")])
# dat.15_45=dat.15_45[col.list]

#subtraction
dat.15_45=dat[1:hc] %>%
  cbind(dat[c("IgG1","TotalG")]) %>%
  cbind(dat.sqrt[c("IgG2ac","IgG2b","IgG3")]) %>%
  cbind(dat.sq2["IgM"])
dat.15_45=dat.15_45[col.list]

#visualize final data
par(mfrow=c(2,3))
for(i in 1:6)
{phenotype<-dat.15_45[,i+hc];hist(phenotype, main=paste0("Kinetics: Day 15 to 45 ",abs[i]),col="honeydew")}
write.csv(dat.15_45,file=paste0("data_bin/kinetics/dat.15-45_",kin.method,"_auc.csv"),row.names=F)



## direction of change
kin.sign=cbind(kin[1:2],sign(kin[3:8]))
kin.method="sign"
dat.7_10=subset(kin.sign,kin.sign$day=="d7_d10")
write.csv(dat.7_10,file=paste0("data_bin/kinetics/dat.7-10_",kin.method,"_auc.csv"),row.names=F)

dat.10_15=subset(kin.sign,kin.sign$day=="d10_d15")
write.csv(dat.7_10,file=paste0("data_bin/kinetics/dat.10-15_",kin.method,"_auc.csv"),row.names=F)

dat.15_45=subset(kin.sign,kin.sign$day=="d15_d45")
write.csv(dat.15_45,file=paste0("data_bin/kinetics/dat.15-45_",kin.method,"_auc.csv"),row.names=F)


rm(dat.inv,dat.inv2,dat.invsqrt,dat.sq2,dat.sqrt,dat.log)
