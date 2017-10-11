#### setup: load libraries, data format, etc ####
setwd("~/Dropbox/Heise/U19-Ab/antibody")

#load libraries 
library(ggplot2)
library(MASS)

#load in wide format data to map
full.dat=read.csv("data_bin/AUC_by_isotype.csv")
data.type="auc"

#set list of anitbody isotypes in the column order they appear in
abs=c("IgG1" ,"IgG2ac", "IgG2b", "IgG3", "IgM", "TotalG")
col.list=colnames(full.dat)

#number of header columns before data
hc=10

#test for zeros
is.zero=function(x){
  ifelse(x<=0,T,F)
}

zero.test=apply(full.dat[11:16],c(1,2),is.zero)
table(zero.test)
full.dat[which(zero.test==T),]


###### visualize data and transform as necessary (individually by timepoint) #####

#### day 7 ####
dat=subset(full.dat,full.dat$day==7)

  # baseline distrubtion
  par(mfrow = c(2,3))
  hist.ab(dat,hc)
  
  #some AUC values are every so slightly below 0.. box cox needs all positive values
  #change the values to be positive non-zero numbers so box-cox will run
  # box.ab(dat,hc=hc,num=6)
  for(i in 1:6)
  {
    phenotype<-dat[,i+hc]
    phenotype[phenotype<=0] <- min(subset(phenotype,phenotype>0),na.rm=T)
    lm.bc=lm(phenotype~RIX,data=dat)
    bc=boxcox(lm.bc)
    print(bc$x[which.max(bc$y)])
  }
  
  #test various transformations
  dat.log<-log.ab(dat,hc)
    colnames(dat.log)<-col.list
    hist.ab(dat.log,hc)
  
  dat.sqrt=trans.ab(dat,0.5,hc)
    colnames(dat.sqrt)<-col.list
    hist.ab(dat.sqrt,hc)
  
  ##D7 Transformations
  #log transform all
  dat.7=dat[1:hc] %>%
    cbind(dat.log[c("IgG1","IgG2ac","IgG2b","IgG3","IgM","TotalG")])
  dat.7=dat.7[col.list]
    
  #visualize final data
  for(i in 1:6)
  {phenotype<-dat.7[,i+hc];hist(phenotype, main=paste0("Day 7 ",abs[i]),col="lavender")}
  write.csv(dat.7,file=paste("data_bin/d7/dat.7.auc.csv",sep=""),row.names=F)
  



#### day 10 ####
dat=subset(full.dat,full.dat$day==10)

  #baseline distribution
  par(mfrow = c(2,3))
  hist.ab(dat,hc)
  
  #box cox
  for(i in 1:6)
  {
    phenotype<-dat[,i+hc]
    phenotype[phenotype<=0] <- min(subset(phenotype,phenotype>0),na.rm=T)
    lm.bc=lm(phenotype~RIX,data=dat)
    bc=boxcox(lm.bc)
    print(bc$x[which.max(bc$y)])
  }
  
  
  #test various transformations
  dat.log<-log.ab(dat,hc)
    colnames(dat.log)<-col.list
    hist.ab(dat.log,hc)
  
  dat.sqrt=trans.ab(dat,0.5,hc)
    colnames(dat.sqrt)<-col.list
    hist.ab(dat.sqrt,hc)
  
  
  #D10 transformations
  #IgG1 = log
  #IgG2ac = sqrt
  #IgG2b = sqrt
  #IgG3 = log
  #IgM = log
  #TotalG = none
  dat.10=dat[1:hc] %>%
    cbind(dat["TotalG"]) %>%  
    cbind(dat.log[c("IgG1","IgG3","IgM")]) %>%
    cbind(dat.sqrt[c("IgG2ac","IgG2b")])
  dat.10=dat.10[col.list]
  
  #visualize final data
  for(i in 1:6)
  {phenotype<-dat.10[,i+hc];hist(phenotype, main=paste0("Day 10 ",abs[i]),col="lavender")}
  write.csv(dat.10,file=paste("data_bin/d10/dat.10.auc.csv",sep=""),row.names=F)

  
  
  
#### day 15 ####
dat=subset(full.dat,full.dat$day==15)

  #baseline distribution
  par(mfrow = c(2,3))
  hist.ab(dat,hc)
  
  for(i in 1:6)
  {
    phenotype<-dat[,i+hc]
    phenotype[phenotype<=0] <- min(subset(phenotype,phenotype>0),na.rm=T)
    lm.bc=lm(phenotype~RIX,data=dat)
    bc=boxcox(lm.bc)
    print(bc$x[which.max(bc$y)])
  }
  
  #test various transformations
  dat.log<-log.ab(dat,hc)
  colnames(dat.log)<-col.list
  hist.ab(dat.log,hc)
  
  dat.sqrt=trans.ab(dat,0.5,hc)
  colnames(dat.sqrt)<-col.list
  hist.ab(dat.sqrt,hc)
  
  
  ##D15 Transformations
  #IgG1 = sqrt
  #IgG2ac = none
  #IgG2b = sqrt
  #IgG3 = sqrt
  #IgM = log
  #TotalG = none
  
  dat.15=dat[1:hc] %>%
    cbind(dat[c("IgG2ac","TotalG")]) %>%
    cbind(dat.log["IgM"]) %>%
    cbind(dat.sqrt[c("IgG1","IgG2b","IgG3")])
  dat.15=dat.15[col.list]
  
  #visualize final data
  for(i in 1:6)
  {phenotype<-dat.15[,i+hc];hist(phenotype, main=paste0("Day 15 ",abs[i]),col="lavender")}
  write.csv(dat.15,paste("data_bin/d15/dat.15.auc.csv",sep=""),row.names=F)         
  
  
#### day 45 ####
dat=subset(full.dat,full.dat$day==45)

  #baseline distribution
  par(mfrow = c(2,3))
  hist.ab(dat,hc)
  
  for(i in 1:6)
  {
    phenotype<-dat[,i+hc]
    phenotype[phenotype<=0] <- min(subset(phenotype,phenotype>0),na.rm=T)
    lm.bc=lm(phenotype~RIX,data=dat)
    bc=boxcox(lm.bc)
    print(bc$x[which.max(bc$y)])
  }
  
  #test various transformations
  dat.log<-log.ab(dat,hc)
  colnames(dat.log)<-col.list
  hist.ab(dat.log,hc)
  
  dat.sqrt=trans.ab(dat,0.5,hc)
  colnames(dat.sqrt)<-col.list
  hist.ab(dat.sqrt,hc)
  
  dat.sq2=trans.ab(dat,2,hc)
  colnames(dat.sq2)<-col.list
  hist.ab(dat.sq2,hc)
  
  dat.invsqrt=trans.ab(dat,-0.5,hc)
  colnames(dat.invsqrt)<-col.list
  hist.ab(dat.invsqrt,hc)
  
  #D45 Transformations
  #IgG1 = sqrt
  #IgG2ac = square
  #IgG2b = none
  #IgG3 = sqrt
  #IgM = 1/sqrt
  #TotalG = square
  dat.45=dat[1:hc] %>%
    cbind(dat["IgG2b"]) %>%
    cbind(dat.invsqrt["IgM"]) %>%
    cbind(dat.sqrt[c("IgG1","IgG3")]) %>%
    cbind(dat.sq2[c("IgG2ac","TotalG")])
  dat.45=dat.45[col.list]
  
  #visualize final data
  for(i in 1:6)
  {phenotype<-dat.45[,i+hc];hist(phenotype, main=paste0("Day 45 ",abs[i]),col="lavender")}
  write.csv(dat.45,file=paste("data_bin/d45/dat.45.auc.csv",sep=""),row.names=F)
  
  alldays.combined=rbind(dat.7,dat.10,dat.15,dat.45)
  write.csv(alldays.combined,file=paste0("data_bin/",data.type,"_transformed_alldays.csv"),row.names=F)
  
# ##### tests for normality ####
# dat=subset(full.dat,full.dat$day==45)
# dat=dat.45
# 
# for(i in 1:6)
# {
#   phenotype<-dat[,i+hc]
#   qqnorm(phenotype, main=abs[i])
# }
# 
# 
# sha.p=NULL
# for(i in 1:6)
# {
#   phenotype<-dat[,i+hc]
#   s.p=shapiro.test(phenotype)$p.value
#   sha.p <- rbind(sha.p, s.p)
# }
# rownames(sha.p)=col.list[(hc+1):ncol(dat)]
# sha.p
  
  
  
  
##### LASTPOS #######
  #### setup: load libraries, data format, etc ####
  setwd("~/Dropbox/Heise/U19-Ab/antibody")
  
  #load in wide format data to map
  full.dat=read.csv("data_bin/last_positive_by_isotype.csv")
  data.type="lastpos"
  
  #set list of anitbody isotypes in the column order they appear in
  abs=c("IgG1" ,  "IgG2ac" , "IgG2b" , "IgG3" ,  "IgM"  ,  "TotalG")
  
  col.list=colnames(full.dat)
  
  #number of header columns before data
  hc=10
  
  # #test for zeros
  # is.zero=function(x){
  #   ifelse(x<=0,T,F)
  # }
  # 
  # zero.test=apply(full.dat[11:16],c(1,2),is.zero)
  # table(zero.test)
  # full.dat[which(zero.test==T),]
  
  
  
  ###### visualize data and transform as necessary (individually by timepoint) #####
  
  #### day 7 ####
  dat=subset(full.dat,full.dat$day==7)
  
  #view distribution
  apply(dat[colnames(dat) %in% abs],2,table)
    
  # baseline distrubtion
  par(mfrow = c(2,3))
  hist.ab(dat,hc)
  
  
  ## not going to normalize binary data - only for totalG ##
  #box cox
  for(i in 1:6)
  {
    phenotype<-dat[,i+hc]
    lm.bc=lm(phenotype~RIX,data=dat)
    bc=boxcox(lm.bc)
    print(bc$x[which.max(bc$y)])
  }

  #test various transformations
  dat.log<-log.ab(dat,hc)
  colnames(dat.log)<-col.list
  hist.ab(dat.log,hc)
  # 
  # dat.inv=trans.ab(dat,-1,hc)
  # colnames(dat.inv)<-col.list
  # hist.ab(dat.inv,hc)
  
  # dat.invsqrt=trans.ab(dat,-0.5,hc)
  # colnames(dat.invsqrt)<-col.list
  # hist.ab(dat.invsqrt,hc)
  
  # dat.invsq2=trans.ab(dat,-2,hc)
  # colnames(dat.invsq2)<-col.list
  # hist.ab(dat.invsq2,hc)
  # 
  ##D7 Transformations

  #binary classification
  dat.bin=dat[c("IgG1","IgG2ac","IgG2b","IgG3","IgM")]
  dat.bin[dat.bin==10] = 0
  dat.bin[dat.bin>=100] = 1
  
  dat.7=dat[1:hc] %>%
    cbind(dat.bin[c("IgG1","IgG2ac","IgG2b","IgG3","IgM")]) %>%
    cbind(dat.log[c("TotalG")]) 
  dat.7=dat.7[col.list]

  
  #visualize final data
  for(i in 1:6)
  {phenotype<-dat.7[,i+hc];hist(phenotype, main=paste0("Day 7 ",abs[i]),col="lavender")}
  write.csv(dat.7,file=paste("data_bin/d7/dat.7.lastpos.csv",sep=""),row.names=F)
  

  
  #### day 10 ####
  dat=subset(full.dat,full.dat$day==10)
  
  #baseline distribution
  par(mfrow = c(2,3))
  hist.ab(dat,hc)
  apply(dat[colnames(dat) %in% abs],2,table)
  
  #box cox
  for(i in 1:6)
  {
    phenotype<-dat[,i+hc]
    phenotype[phenotype<=0] <- min(subset(phenotype,phenotype>0),na.rm=T)
    lm.bc=lm(phenotype~RIX,data=dat)
    bc=boxcox(lm.bc)
    print(bc$x[which.max(bc$y)])
  }
  
  #test various transformations
  dat.log<-log.ab(dat,hc)
  colnames(dat.log)<-col.list
  hist.ab(dat.log,hc)

  #D10 transformations
  #log all but IgG1
  dat.bin=dat["IgG1"]
  dat.bin[dat.bin==10] = 0
  dat.bin[dat.bin>=100] = 1
  
  dat.10=dat[1:hc] %>%
    cbind(dat.bin) %>%
    cbind(dat.log[c("IgG2ac","IgG2b","IgG3","IgM","TotalG")])
  dat.10=dat.10[col.list]
  
  #visualize final data
  for(i in 1:6)
  {phenotype<-dat.10[,i+hc];hist(phenotype, main=paste0("Day 10 ",abs[i]),col="lavender")}
  write.csv(dat.10,file=paste("data_bin/d10/dat.10.lastpos.csv",sep=""),row.names=F)
  
  
  
  
  #### day 15 ####
  dat=subset(full.dat,full.dat$day==15)
  
  #baseline distribution
  par(mfrow = c(2,3))
  hist.ab(dat,hc)
  apply(dat[colnames(dat) %in% abs],2,table)
  
  for(i in 1:6)
  {
    phenotype<-dat[,i+hc]
    phenotype[phenotype<=0] <- min(subset(phenotype,phenotype>0),na.rm=T)
    lm.bc=lm(phenotype~RIX,data=dat)
    bc=boxcox(lm.bc)
    print(bc$x[which.max(bc$y)])
  }
  
  #test various transformations
  dat.log<-log.ab(dat,hc)
  colnames(dat.log)<-col.list
  hist.ab(dat.log,hc)
  
  
  ##D15 Transformations
  #binary IgM, log the rest
  dat.bin=dat["IgM"]
  dat.bin[dat.bin==10] = 0
  dat.bin[dat.bin>=100] = 1
  
  dat.15=dat[1:hc] %>%
    cbind(dat.bin) %>%
    cbind(dat.log[c("IgG1","IgG2ac","IgG2b","IgG3","TotalG")])
  dat.15=dat.15[col.list]
  
  #visualize final data
  for(i in 1:6)
  {phenotype<-dat.15[,i+hc];hist(phenotype, main=paste0("Day 15 ",abs[i]),col="lavender")}
  write.csv(dat.15,paste("data_bin/d15/dat.15.lastpos.csv",sep=""),row.names=F)         
  
  
  #### day 45 ####
  dat=subset(full.dat,full.dat$day==45)
  
  #baseline distribution
  par(mfrow = c(2,3))
  hist.ab(dat,hc)
  apply(dat[colnames(dat) %in% abs],2,table)
  
  for(i in 1:6)
  {
    phenotype<-dat[,i+hc]
    phenotype[phenotype<=0] <- min(subset(phenotype,phenotype>0),na.rm=T)
    lm.bc=lm(phenotype~RIX,data=dat)
    bc=boxcox(lm.bc)
    print(bc$x[which.max(bc$y)])
  }
  
  #test various transformations
  dat.log<-log.ab(dat,hc)
  colnames(dat.log)<-col.list
  hist.ab(dat.log,hc)
  
  dat.sqrt=trans.ab(dat,0.5,hc)
  colnames(dat.sqrt)<-col.list
  hist.ab(dat.sqrt,hc)
  
  dat.sq2=trans.ab(dat,2,hc)
  colnames(dat.sq2)<-col.list
  hist.ab(dat.sq2,hc)
  
  dat.invsqrt=trans.ab(dat,-0.5,hc)
  colnames(dat.invsqrt)<-col.list
  hist.ab(dat.invsqrt,hc)
  
  #D45 Transformations
  #IgG1 = log
  #IgG2ac = none
  #IgG2b = sqrt
  #IgG3 = log
  #IgM = binary
  #TotalG = square
  
  dat.bin=dat["IgM"]
  dat.bin[dat.bin==10] = 0
  dat.bin[dat.bin>=100] = 1
  
  dat.45=dat[1:hc] %>%
    cbind(dat["IgG2ac"]) %>%
    cbind(dat.bin["IgM"]) %>%
    cbind(dat.sqrt[c("IgG2b")]) %>%
    cbind(dat.log[c("IgG1","IgG3")]) %>%
    cbind(dat.sq2["TotalG"])
  dat.45=dat.45[col.list]
  
  #visualize final data
  for(i in 1:6)
  {phenotype<-dat.45[,i+hc];hist(phenotype, main=paste0("Day 45 ",abs[i]),col="lavender")}
  write.csv(dat.45,file=paste("data_bin/d45/dat.45.lastpos.csv",sep=""),row.names=F)
  
  alldays.combined=rbind(dat.7,dat.10,dat.15,dat.45)
  write.csv(alldays.combined,file=paste0("data_bin/",data.type,"_transformed_alldays.csv"),row.names=F)
  
##########
  
  ##### HALFMAX #######
  #### setup: load libraries, data format, etc ####
  setwd("~/Dropbox/Heise/U19-Ab/antibody")
  data.type="halfmax"
  
  #load in wide format data to map
  full.dat=read.csv("data_bin/halfmax_by_isotype.csv")
  
  #set list of anitbody isotypes in the column order they appear in
  abs=c("IgG1" ,  "IgG2ac" , "IgG2b" , "IgG3" ,  "IgM"  ,  "TotalG")
  
  col.list=colnames(full.dat)
  
  #number of header columns before data
  hc=10
  
  # #test for zeros
  # is.zero=function(x){
  #   ifelse(x<=0,T,F)
  # }
  # 
  # zero.test=apply(full.dat[11:16],c(1,2),is.zero)
  # table(zero.test)
  # full.dat[which(zero.test==T),]
  
  
  
  ###### visualize data and transform as necessary (individually by timepoint) #####
  
  #### day 7 ####
  dat=subset(full.dat,full.dat$day==7)
  
  #view distribution
  apply(dat[colnames(dat) %in% abs],2,table)
  
  # baseline distrubtion
  par(mfrow = c(2,3))
  hist.ab(dat,hc)
  
  #binary classification
  dat.bin=dat[c("IgG1","IgG2ac","IgG2b","IgG3","IgM","TotalG")]
  dat.bin[dat.bin==10] = 0
  dat.bin[dat.bin>=100] = 1
  
  dat.7=dat[1:hc] %>%
    cbind(dat.bin[c("IgG1","IgG2ac","IgG2b","IgG3","IgM","TotalG")])
  dat.7=dat.7[col.list]
  
  
  #visualize final data
  for(i in 1:6)
  {phenotype<-dat.7[,i+hc];hist(phenotype, main=paste0("Day 7 ",abs[i]),col="lavender")}
  write.csv(dat.7,file=paste0("data_bin/d7/dat.7.",data.type,".csv"),row.names=F)
  
  
  #### day 10 ####
  dat=subset(full.dat,full.dat$day==10)
  
  #baseline distribution
  par(mfrow = c(2,3))
  hist.ab(dat,hc)
  apply(dat[colnames(dat) %in% abs],2,table)
  
  #box cox
  for(i in 1:6)
  {
    phenotype<-dat[,i+hc]
    phenotype[phenotype<=0] <- min(subset(phenotype,phenotype>0),na.rm=T)
    lm.bc=lm(phenotype~RIX,data=dat)
    bc=boxcox(lm.bc)
    print(bc$x[which.max(bc$y)])
  }
  
  #test various transformations
  dat.log<-log.ab(dat,hc)
  colnames(dat.log)<-col.list
  hist.ab(dat.log,hc)
  
  #D10 transformations
  #log all IgG2ac and TotalG
  #binary the rest
  dat.bin=dat[c("IgG1","IgG2b","IgG3","IgM")]
  dat.bin[dat.bin==10] = 0
  dat.bin[dat.bin>=100] = 1
  
  dat.10=dat[1:hc] %>%
    cbind(dat.bin[c("IgG1","IgG2b","IgG3","IgM")]) %>%
    cbind(dat.log[c("IgG2ac","TotalG")])
  dat.10=dat.10[col.list]
  
  #visualize final data
  for(i in 1:6)
  {phenotype<-dat.10[,i+hc];hist(phenotype, main=paste0("Day 10 ",abs[i]),col="lavender")}
  write.csv(dat.10,file=paste0("data_bin/d10/dat.10.",data.type,".csv"),row.names=F)
  
  
  #### day 15 ####
  dat=subset(full.dat,full.dat$day==15)
  
  #baseline distribution
  par(mfrow = c(2,3))
  hist.ab(dat,hc)
  apply(dat[colnames(dat) %in% abs],2,table)
  
  for(i in 1:6)
  {
    phenotype<-dat[,i+hc]
    phenotype[phenotype<=0] <- min(subset(phenotype,phenotype>0),na.rm=T)
    lm.bc=lm(phenotype~RIX,data=dat)
    bc=boxcox(lm.bc)
    print(bc$x[which.max(bc$y)])
  }
  
  #test various transformations
  dat.log<-log.ab(dat,hc)
  colnames(dat.log)<-col.list
  hist.ab(dat.log,hc)
  
  
  ##D15 Transformations
  #binary IgM, log the rest
  dat.bin=dat[c("IgG1","IgG2b","IgG3","IgM")]
  dat.bin[dat.bin==10] = 0
  dat.bin[dat.bin>=100] = 1
  
  dat.15=dat[1:hc] %>%
    cbind(dat.bin) %>%
    cbind(dat.log[c("IgG2ac","TotalG")])
  dat.15=dat.15[col.list]
  
  #visualize final data
  for(i in 1:6)
  {phenotype<-dat.15[,i+hc];hist(phenotype, main=paste0("Day 15 ",abs[i]),col="lavender")}
  write.csv(dat.15,paste0("data_bin/d15/dat.15.",data.type,".csv"),row.names=F)         
  
  #### day 45 ####
  dat=subset(full.dat,full.dat$day==45)
  
  #baseline distribution
  par(mfrow = c(2,3))
  hist.ab(dat,hc)
  apply(dat[colnames(dat) %in% abs],2,table)
  
  for(i in 1:6)
  {
    phenotype<-dat[,i+hc]
    phenotype[phenotype<=0] <- min(subset(phenotype,phenotype>0),na.rm=T)
    lm.bc=lm(phenotype~RIX,data=dat)
    bc=boxcox(lm.bc)
    print(bc$x[which.max(bc$y)])
  }
  
  #test various transformations
  dat.log<-log.ab(dat,hc)
  colnames(dat.log)<-col.list
  hist.ab(dat.log,hc)
  
  dat.sqrt=trans.ab(dat,0.5,hc)
  colnames(dat.sqrt)<-col.list
  hist.ab(dat.sqrt,hc)
  
  #D45 Transformations
  #IgG1 = log
  #IgG2ac = none
  #IgG2b = sqrt
  #IgG3 = log
  #IgM = binary
  #TotalG = square
  
  dat.bin=dat[c("IgG1","IgG3","IgM")]
  dat.bin[dat.bin==10] = 0
  dat.bin[dat.bin>=100] = 1
  
  dat.45=dat[1:hc] %>%
    cbind(dat.bin[c("IgG1","IgG3","IgM")]) %>%
    cbind(dat.log[c("IgG2ac","IgG2b","TotalG")])
  dat.45=dat.45[col.list]
  
  #visualize final data
  for(i in 1:6)
  {phenotype<-dat.45[,i+hc];hist(phenotype, main=paste0("Day 45 ",abs[i]),col="lavender")}
  write.csv(dat.45,file=paste0("data_bin/d45/dat.45.",data.type,".csv"),row.names=F)
  
  alldays.combined=rbind(dat.7,dat.10,dat.15,dat.45)
  write.csv(alldays.combined,file=paste0("data_bin/",data.type,"_transformed_alldays.csv"),row.names=F)
  
  
  # ##### tests for normality ####
  # dat=subset(full.dat,full.dat$day==45)
  # dat=dat.45
  # 
  # for(i in 1:6)
  # {
  #   phenotype<-dat[,i+hc]
  #   qqnorm(phenotype, main=abs[i])
  # }
  # 
  # 
  # sha.p=NULL
  # for(i in 1:6)
  # {
  #   phenotype<-dat[,i+hc]
  #   s.p=shapiro.test(phenotype)$p.value
  #   sha.p <- rbind(sha.p, s.p)
  # }
  # rownames(sha.p)=col.list[(hc+1):ncol(dat)]
  # sha.p
  
  
  
  
  