#### SETUP: load libraries, data format, etc ####
setwd("~/Dropbox/Heise/U19-Ab/sars_ab/")

#load libraries 
library(ggplot2)
library(MASS)

#select dataset (lastpos or halfmax)
# dataset="lastpos"
dataset="halfmax"

#select antigen (N or S)
anti="S"

#load in data to map in wide format
full.dat=read.csv(paste0("data_bin/SARS_",dataset,"_by_isotype.csv"))

#set number of ID header columns before data
hc=6

#set proper order of days as factor levels
full.dat$day = full.dat$day %>% relevel("DRECHAL") %>% relevel("D29") %>% relevel("D15") %>% relevel("D10") %>% relevel("D7")

#capture list of column and antibody isotype names as they appear in data
col.list=colnames(full.dat)
abs=colnames(full.dat[(hc+1):length(full.dat)])




##### visualize data and transform as necessary (individually by timepoint#) #####
dat.anti=subset(full.dat,full.dat$antigen==anti)

#### day ####
day="D7"

dat=subset(dat.anti,dat.anti$day==day)

# graph to view data distribution without transformation #
par(mfrow = c(2,3))
for(i in 1:6)
{
  phenotype<-dat[,i+hc]
  hist(phenotype, main=abs[i])
}

#run box cox analysis to get an idea of best transformations
#box cox can't deal with 0s or negatives so if you have them in your data frame you need to adjust
for(i in 1:6)
{
  phenotype<-dat[,i+hc]
  lm.bc=lm(phenotype~RIX,data=dat)
  bc=boxcox(lm.bc)
  print(bc$x[which.max(bc$y)])
}

#see what they look like log transformed - more normal?
dat.log<-log.ab(dat,hc)
  colnames(dat.log)<-col.list

    for(i in 1:6)    {
      phenotype<-dat.log[,i+hc]
      hist(phenotype, main=abs[i])
    }

#see what they look like sqrt transformed - more normal?
dat.sqrt=trans.ab(dat,0.5,hc)
  colnames(dat.sqrt)<-col.list
    
    for(i in 1:6)    {
      phenotype<-dat.sqrt[,i+hc]
      hist(phenotype, main=abs[i])
     }

#see what they look like inversed - more normal?
dat.inv=trans.ab(dat,-1,hc)
  colnames(dat.inv)<-col.list
    
    for(i in 1:6) {
      phenotype<-dat.inv[,i+hc]
      hist(phenotype, main=abs[i])
     }    
    
#see what they look like sq.inversed - more normal?
dat.inv2<-trans.ab(dat,-2,hc)
  colnames(dat.inv2)<-col.list
    
    for(i in 1:6)  {
      phenotype<-dat.inv2[,i+hc]
      hist(phenotype, main=abs[i])
     }    

  
  
#### D7 Transformations ####
dat.7=dat[1:hc] %>%
  cbind(dat["IgG1"]) %>%
  cbind(dat.log[c("IgG2ac","IgM","TotalG")]) %>%
  cbind(dat.inv["IgG2b"]) %>%
  cbind(dat.inv2["IgG3"])
dat.7=dat.7[col.list]

#visualize and export final data
for(i in 1:6)
{
  phenotype<-dat.7[,i+hc]
  hist(phenotype, main=abs[i],col="lavender")
}
write.csv(dat.7,file=paste0("data_bin/SARS_D7_",anti,"_",dataset,".csv"),row.names=F)




#### D10 Transformations ####
#IgG1 = log
#IgG2ac = sqrt
#IgG2b = long
#IgG3 = log
#IgM = log
#TotalG = sqrt
dat.10=dat[1:hc] %>%
  cbind(dat.log[c("IgG1","IgG3","IgM")]) %>%
  cbind(dat.sqrt[c("IgG2ac","IgG2b","TotalG")])
dat.10=dat.10[col.list]

#visualize final data
for(i in 1:6)
{
  phenotype<-dat.10[,i+hc]
  hist(phenotype, main=abs[i],col="lavender")
}

write.csv(dat.10,file=paste0("data_bin/SARS_D10_",anti,"_",dataset,".csv"),row.names=F)





#### D15 Transformations ####
#IgG1 = sqrt
#IgG2ac = sqrt
#IgG2b = sqrt
#IgG3 = log
#IgM = log
#TotalG = none

dat.15=dat[1:hc] %>%
  cbind(dat[c("TotalG")]) %>%
  cbind(dat.log[c("IgG3","IgM")]) %>%
  cbind(dat.sqrt[c("IgG1","IgG2ac","IgG2b")])
dat.15=dat.15[col.list]

#visualize final data
for(i in 1:6)
{
  phenotype<-dat.15[,i+hc]
  hist(phenotype, main=abs[i],col="lavender")
}

write.csv(dat.15,file=paste0("data_bin/SARS_D15_",anti,"_",dataset,".csv"),row.names=F)





#### D29 Transformations ####
#IgG1 = sqrt
#IgG2ac = sqrt
#IgG2b = sqrt
#IgG3 = log
#IgM = log
#TotalG = none

dat.29=dat[1:hc] %>%
  cbind(dat[c("TotalG")]) %>%
  cbind(dat.log[c("IgG3","IgM")]) %>%
  cbind(dat.sqrt[c("IgG1","IgG2ac","IgG2b")])
dat.29=dat.29[col.list]

#visualize final data
for(i in 1:6)
{
  phenotype<-dat.29[,i+hc]
  hist(phenotype, main=abs[i],col="lavender")
}

write.csv(dat.29,file=paste0("data_bin/SARS_D29_",anti,"_",dataset,".csv"),row.names=F)





#### Rechallenge ####
##Drechal Transformations
#IgG1 = sqrt
#IgG2ac = sqrt
#IgG2b = sqrt
#IgG3 = log
#IgM = log
#TotalG = none

dat.rechal=dat[1:hc] %>%
    cbind(dat[c("TotalG")]) %>%
    cbind(dat.log[c("IgG3","IgM")]) %>%
    cbind(dat.sqrt[c("IgG1","IgG2ac","IgG2b")])
dat.rechal=dat.rechal[col.list]

write.csv(dat.rechal,file=paste0("data_bin/SARS_DRECHAL_",anti,"_",dataset,".csv"),row.names=F)
          