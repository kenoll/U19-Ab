#### SETUP ####

# load required packages
library(dplyr)
library(reshape2)
library(tidyr)
library(flux)
library("data.table")

#set working directory
setwd("~/Dropbox/Heise/U19-Ab/sars_ab/")



#### INPUT DATA ####

#load in data file in wide plate format with first column as dilution value
#fread function in data.table package reads in faster than read.csv
data1000 <- as.data.frame(fread("data_bin/SARS_1_1000.csv",stringsAsFactors=FALSE, header=TRUE,sep=","))
data100 <- as.data.frame(fread("data_bin/SARS_1_100.csv",stringsAsFactors=FALSE, header=TRUE,sep=","))

#convert wide plate format to long table
data100=plates.to.table(data100)
data1000=plates.to.table(data1000)

#add in empty columns for missing dilutions to 100 and 1000 data frames so they are the same across both
#bind data frames
data100$X_1000000=NA
data100$X_3000000=NA
data1000$X_100=NA
data1000$X_300=NA
data1=rbind(data100,data1000)

#provide information for data table structure
ndil=10 #number of different dilutions
hc=7    #number of 'header columns' (sample ID and isotype information before absorbance values)

#create table index of dilution factors used and corresponding indeces
#the first value "1" is an arbitrary placeholder used later for a value below the limit of detection
dilutions=c("1","100","300","1000","3000","10000","30000","100000","300000","1000000","3000000")
dilutions.index=c(0,1,2,3,4,5,6,7,8,9,10)
dilframe=data.frame(dilutions.index,dilutions)

#separate sample information into multiple columns
data1=separate(data1,ID,into=c("RIX","ID","day","isotype"),sep="_",extra="merge")
data1=separate(data1,isotype,into=c("isotype","virus","antigen"),sep="_")

#remove B6 control sera into a separate data frame
controls=data1[data1$ID %in% "B6",]
data1=anti_join(data1,controls)

#convert to new CC names
data1$alias=data1$RIX
data1=data1[c(1,length(data1),2:(length(data1)-1))]
data1=alias.to.line(data1)

#convert absorbance values to numeric
data1[(hc+1):length(data1)] = sapply(data1[(hc+1):length(data1)], as.character)
data1[(hc+1):length(data1)] = sapply(data1[(hc+1):length(data1)], as.numeric)

#remove duplicates (repeated samples) until clarified
dups=which(duplicated(data1[c("RIX","ID","day","isotype","antigen")]) |
                duplicated(data1[c("RIX","ID","day","isotype","antigen")], fromLast = TRUE))
dups=data1[dups,]
dups=dups[order(dups$RIX,dups$ID,dups$day,dups$isotype,dups$antigen),]
data1=anti_join(data1,dups)




#### HALF MAX (absorbance >1.75) ####

#function takes anything under the threshold and sets it to 5 (arbitrary high number)
#then uses which.min to pick out the index of the lowest value in that column
#note: outputs INDEX of row, not actual log dilution value
dilfind=function(x){
  x[x<1.75] = 5
  which.min(x)
  }

#apply dilfind function over each column of the data frame
dil=NULL
dil=apply(data1[(ncol(data1)-ndil+1):(ncol(data1))],1,dilfind)

#if the row with the lowest dilution has a value less than threshold, replace that row index with 0
#using "0" instead of NA to distinguish between sample not run
for (i in 1:nrow(data1)){
  if(data1[i,(dil[[i]]+hc)]<1.75)
  {dil[[i]]=0}}

#create new dataframe with halfmax values
data1[(ncol(data1)+1)]=dil
colnames(data1)[(ncol(data1))]="halfmax"
halfmax=data1[c(1:hc,ncol(data1))]

#drop halfmax from absorbance data
data1=data1[-ncol(data1)]

#convert halfmax from column index to dilution factor
colnames(dilframe)[1]=c("halfmax")
halfmax=merge(halfmax,dilframe,all=T)
halfmax=halfmax[-1]
colnames(halfmax)[length(halfmax)]=c("halfmax")
halfmax = halfmax[order(halfmax$RIX),]

#convert halfmax data to numeric
halfmax$halfmax= as.character(halfmax$halfmax)
halfmax$halfmax= as.numeric(halfmax$halfmax)

#cast back into wide format based on isotype
wide.halfmax=dcast(halfmax, RIX + alias + ID + day + virus + antigen ~ isotype,mean,value.var='halfmax')
write.csv(wide.halfmax,"data_bin/SARS_halfmax_by_isotype.csv",row.names=F)



#### LAST POSITIVE (absorbance > 0.2) ####

#function takes anything under the threshold and sets it to 5 (arbitrary high number)
#then uses which.min to pick out the index of the lowest value in that column
#note: outputs INDEX of row, not actual log dilution value
dilfind=function(x){x[x<0.2] = 5;which.min(x)}

#apply over each column of the data frame
dil=NULL
dil=apply(data1[(ncol(data1)-ndil+1):(ncol(data1))],1,dilfind)

#if the row with the lowest dilution has a value less than threshold, replace that row index with 0
#using "0" instead of NA to distinguish between sample not run
for (i in 1:nrow(data1)){
  if(data1[i,(dil[[i]]+hc)]<0.2)
  {dil[[i]]=0}}

#create new data frame with last positive data
data1[(ncol(data1)+1)]=dil
colnames(data1)[ncol(data1)]="last_positive"
lastpos=data1[c(1:hc,ncol(data1))]

#drop last positive from absorbance data
data1=data1[-ncol(data1)]

#convert last positive value from column index to dilution factor
colnames(dilframe)[1]=c("last_positive")
lastpos=merge(lastpos,dilframe,all=T)
lastpos=lastpos[-1]
colnames(lastpos)[length(lastpos)]=c("last_positive")

#convert last positive value to numeric
lastpos$last_positive= as.character(lastpos$last_positive)
lastpos$last_positive= as.numeric(lastpos$last_positive)

#cast back into wide format based on isotype
wide.lastpos=dcast(lastpos, RIX + alias + ID + day + virus + antigen ~ isotype,mean,value.var='last_positive')
write.csv(wide.lastpos,paste0("data_bin/SARS_lastpos_by_isotype.csv"),row.names=F)



#### THE GRANDE FINALE ####

#merge long form last positive and half max data
hm.lp=merge(halfmax,lastpos)
write.csv(hm.lp,"data_bin/SARS_halfmax_lastpos.csv",row.names=F)

#clean environment
rm(data1,data100,data1000,ndil,dil,dilframe,dilutions,dilutions.index,dups)