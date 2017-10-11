#### loads ####
library(reshape2)
library(tidyr)
library(flux)
library(dplyr)
library(doBy)
library(data.table)

#### getting started ####
#load file and make first column (dilutions) the row name
#the loaded data file has the first colum as the dilution vales
#the bottom row is the dates on which the ELISAs were run
#samples with wonky curves have been removed manually already (as per QC graphs at the end of the file)
setwd("~/Dropbox/Heise/U19-Ab/antibody")

# data1=read.csv("data_bin/raw_platereader_2017.csv")
data1<-as.data.frame(fread("data_bin/raw_platereader_2017.csv",stringsAsFactors=FALSE, header=TRUE,sep=","))
  data1=plates.to.table(data1,date=T)
  data1=separate(data1,ID,into=c("RIX","ID","day","isotype"),sep="_",extra="merge")
  data1=separate(data1,isotype,into=c("isotype","antigen_virus","antigen"),sep="_")
    #will get a warning - ignore! :)
  data1$antigen_virus="Influenza"
  data1$antigen="HA"

#dilutions
  dilutions=c("10","100","300","1000","3000","10000","30000","100000","300000")
  dilutions.index=c(0,1,2,3,4,5,6,7,8)
  dilframe=data.frame(dilutions.index,dilutions)
  
  
#### split RIX ####
  data1=split.rix(data1)
  data1$dam=gsub("X","",data1$dam)

#set header columns
  hc=9

  ### cleanup data
  #remove 16034x13067_45_0 from data bc it's a mock and it has a duplicate ID number that's messing stuff up
  data1=data1[-which(data1$RIX=="X16034x13067" & data1$ID=="45" & data1$day=="0"),]

  ## get rid of or fix day for 8016x8004_30 (alan wrote in as 0, should be 15 - there already is a sample labeled 15)
  # data1[which(data1$RIX=="X8016x8004" & data1$ID=="30" & data1$day=="0"),]$day="15"
  data1=data1[-which(data1$RIX=="X8016x8004" & data1$ID=="30" & data1$day=="0"),]
  # 
  # # #IgG2ab are replicates of IgG2ac done on the same day, remove them and don't lose anything
  # data1=data1[-which(data1$isotype=="IgG2ab"),]
  
  #IgG2ab = IgG2ac
  data1[which(data1$isotype=="IgG2ab"),]$isotype="IgG2ac"
  
  #fix misentered alias names
  data1[which(data1$dam=="13067" & data1$sire=="8306"),]$sire="5306"
  data1[which(data1$dam=="16411" & data1$sire=="8005"),]$dam="16441"
  data1[which(data1$dam=="5035" & data1$sire=="16875"),]$sire="16785"
  data1[which(data1$dam=="5306" & data1$sire=="5436"),]$sire="5346"
  data1[which(data1$dam=="8002" & data1$sire=="8032"),]$sire="3032"
  
  #take out controls/negatives/etc
  data1=subset(data1,data1$RIX!="negxneg" & data1$RIX!="B6xB6" & data1$RIX!="B6.xB6." & data1$RIX!="X129x129" & data1$RIX!="antixanti")
  
  #remove samples that have NA across all dilutions (failed run)
  ind <- apply(data1[(hc+1):(hc+8)], 1, function(x) all(is.na(x)))
    data1 <- data1[ !ind, ]
    
  #get rid of the d0 mocks so can bind by day
  mocks=data1[which(data1$day=="0"),]
  data1=data1[-which(data1$day=="0"),]
  
  #fix day errors
  data1["X8018x3154_31_102_TotalG",]$day=7
  data1[which(data1$RIX=="X16441x8005" & data1$ID=="76"),]$day=45
  
  #fix sample ID errors
  data1[which(data1$RIX=="X8033x5346" & data1$ID=="55.28." & data1$day=="7"),]$ID=55
  data1[which(data1$RIX=="X8033x5346" & data1$ID=="56.31." & data1$day=="7"),]$ID=56
  
  #examine duplicates
  data1$isotype=gsub("\\.[[:digit:]]","",data1$isotype)
  dups=data1[which(duplicated(data1[c("RIX","ID","day","isotype")]) | 
                      duplicated(data1[c("RIX","ID","day","isotype")], fromLast = TRUE)),]
  dups=dups[order(dups$RIX,dups$ID,dups$isotype,dups$day),]
  
  #remove duplicates from data frame
  data1=anti_join(data1,dups)

        # #average dups if they have the same isotype and were done on the same day (replicate runs)
        # dups=summaryBy(X_2+X_2.5+X_3+X_3.5+X_4+X_4.5+X_5+X_5.5 ~ ., data=dups, FUN=mean, na.rm=T)
        # colnames(dups)=gsub(".mean","",colnames(dups))
        # fixed=dups[which(!duplicated(dups[c("RIX","ID","day","isotype")]) & 
        #                    !duplicated(dups[c("RIX","ID","day","isotype")], fromLast = TRUE)),]
        # 
        # #merge same day replicates with the original data frame and take them away from duplicate df
        # data1=rbind(data1,fixed)
        # dups=anti_join(dups,fixed)
        # 
        # dups$assay_date=droplevels(dups$assay_date)
        # dups=dups[order(dups$RIX,dups$ID,dups$isotype,dups$day),]

  #compare variation between duplicates and merge ones that are close
  dups.n=aggregate(X_2 ~ dam+sire+RIX+ID+day+isotype+antigen_virus+antigen,data=dups,FUN=length)
    colnames(dups.n)[length(dups.n)]="n"
  dups.sd=aggregate(X_2 ~ dam+sire+RIX+ID+day+isotype+antigen_virus+antigen,data=dups,FUN=sd)
    colnames(dups.sd)[length(dups.sd)]="sd"
  dups.var = dups %>% merge(dups.sd) %>% merge(dups.n)
  
  sd.hist=hist(dups.var$sd)
    sd.hist=data.frame(sd.hist$breaks[1:(length(sd.hist$breaks)-1)],sd.hist$counts)
    colnames(sd.hist)=c("sd bin","n")
  
  low.sd=dups[which(dups.var$sd<0.4),]
  low.sd=summaryBy(X_2+X_2.5+X_3+X_3.5+X_4+X_4.5+X_5+X_5.5 ~ 
                     dam+sire+RIX+ID+day+isotype+antigen_virus+antigen, data=low.sd, FUN=mean, na.rm=T)
  colnames(low.sd)=gsub(".mean","",colnames(low.sd))
  low.sd$assay_date="merged"
  
  #return merged close duplicates to data frame
  data1=rbind(data1,low.sd)
  
  # #max absorbance variation correction?
  # od=data1[9:10]
  # od=summaryBy(X_2~assay_date,data = od,FUN = max,na.rm=T)
  # od=subset(od,od$X_2.max>3)
  
  
  # test=dups[which(duplicated(dups[c("RIX","ID","day","isotype")]) | 
                     # duplicated(dups[c("RIX","ID","day","isotype")], fromLast = TRUE)),]
  
#switch to newer CC nomenclature
data1=alias.to.line2(data1)
  
# ### infection status ###
#infection_status.R
  
infections=read.csv("data_bin/infection_status.csv")

# data1=merge(data1,infections[c(1:2,4:6)],all.x=T) #ignore "day" in infection status bc Alan changed a lot of mock days to 0
data1=merge(data1,infections,all.x=T)
data1=data1[c(1:hc,(length(data1)-3):length(data1),(hc+1):(length(data1)-4))]
hc=hc+4
data1[which(data1$day==29),]$virus="Influenza"

#### QC ####
# 
# ##pull out rows w/ no attached infection status to further examine - samples with  no ID #s etc
  mini.dat=data1[c("RIX","dam","sire","ID","day","virus")]
  mini.dat=unique(mini.dat)
  novirus=mini.dat[which(is.na(mini.dat$virus)),]
  novirus=line.to.alias(novirus)

#   mini.dat=data1[c("RIX_ID","RIX","dam","sire","ID","day")]
#   mini.dat=unique(mini.dat)
#   mini.dat$RIX=as.factor(mini.dat$RIX)
#   mini.dat[c("ID","day")]=sapply(mini.dat[c("ID","day")],as.integer)
#   mini.inf=infections[which(infections$day %in% c(7,10,15,29,45 & infections$virus=="Influenza")),]
#   mini.inf=mini.inf[c("RIX_ID","RIX","ID","day")]
#   
#   novirus=anti_join(mini.dat,mini.inf)
#     novirus=subset(novirus,novirus$day!=29)
#   noAb=anti_join(mini.inf,mini.dat)
#   
#   rm(mini.dat,mini.inf)

  
#remove any newly ID'd mocks (not removed as d0) or samples w/o infection status
mocks=subset(data1,data1$virus=="Mock")  
data1=subset(data1,data1$virus=="Influenza")  
  
#make sure all are numeric and 
data1[(length(data1)-7):length(data1)] = sapply(data1[(length(data1)-7):length(data1)], as.character)
data1[(length(data1)-7):length(data1)] = sapply(data1[(length(data1)-7):length(data1)], as.numeric)

#set maximal OD to 3.5 to equalize
data.ods=data1[(length(data1)-7):length(data1)]
data.ods[data.ods>3.5] = 3.504
data1=cbind(data1[1:hc],data.ods)

#subtract background
data1[(length(data1)-7):length(data1)] = data1[(length(data1)-7):length(data1)]-0.04 #generic BG subtraction
data1[, (length(data1)-7):length(data1)][data1[, (length(data1)-7):length(data1)] < 0] = 0 #make negatives zeroes for AUC


#### Calculate AUC ####
#make vector for dilutions to do calculations with later
dilutions=c(1/100,1/300,1/1000,1/3000,1/10000,1/30000,1/100000,1/300000)
dilutions=abs(log10(dilutions))

#make a new row in the df to populate with AUC data
aucdata=data1
aucdata$AUC=NA

#calculate AUC and add to data frame in last row
late=aucdata[which(aucdata$day %in% c("10","15","29","45")),]
late=late[complete.cases(late[(length(aucdata)-8):(length(aucdata)-1)]),]
for(i in 1:nrow(late))
{
  late[i,ncol(late)]<-auc(dilutions, late[i,(ncol(late)-8):(ncol(late)-1)])
}  

seven=aucdata[which(aucdata$day %in% c("7")),]
seven=seven[complete.cases(seven[1:(ncol(seven)-5)]),]
for(i in 1:nrow(seven))
{
  seven[i,ncol(seven)]<-auc(dilutions[1:6], seven[i,(ncol(seven)-8):(ncol(seven)-3)])
}  

aucdata=rbind(late,seven)
rm(late,seven)

#create new DF with sample IDs and AUC data only and export
aucdata=aucdata[c(1:hc,length(aucdata))]
write.csv(aucdata,"data_bin/AUC_data.csv",row.names=F)

# #pull out duplicates to examine for QC
# dup.auc=which(duplicated(aucdata[c("RIX","ID","day","isotype")]) | 
#                 duplicated(aucdata[c("RIX","ID","day","isotype")], fromLast = TRUE))
# aucdata.dup=aucdata[dup.auc,]
# aucdata.dup=aucdata.dup[order(aucdata.dup$RIX,aucdata.dup$ID,aucdata.dup$day,aucdata.dup$isotype),]
# 
# # write.csv(aucdata.dup,"AUC_duplicated.csv",row.names=F)

# cast back into wide format based on isotype
# will average out any duplicate values (e.g. same sample run different days)

wideauc=dcast(aucdata, RIX + dam + sire + ID + day + virus + antigen_virus + antigen + cohort + assay_date ~ isotype,mean,value.var='AUC')
write.csv(wideauc,"data_bin/AUC_by_isotype.csv",row.names=F)


#### HALF MAX (absorbance >1.75) ####

#function takes anything under the threshold and sets it to 5 (arbitrary high number)
#then uses which.min to pick out the index of the lowest value in that column
#note: outputs INDEX of row, not actual log dilution value
dilfind=function(x){
  x[x<1.75] = 5
  unname(which.min(x))
}

#apply over each column of the data frame
dil=NULL
dil=apply(data1[,(hc+1):(hc+8)],1,dilfind)
dil=as.numeric(dil)

#if the row with the lowest dilution has a value less than threshold
#replace that row index with 0
#using "0" instead of NA to distinguish between sample not run
for (i in 1:length(dil)){
  if(data1[i,(dil[[i]]+hc)]<1.75) {
    dil[[i]]=0
  }
  }

#create new halfmax data frame
halfmax=cbind(data1[1:hc],dil)
colnames(halfmax)[(ncol(halfmax))]="halfmax"

#convert halfmax value from column index to dilution factor
colnames(dilframe)[1]=c("halfmax")
halfmax=merge(halfmax,dilframe,all=T)
halfmax=halfmax[-1]
colnames(halfmax)[length(halfmax)]=c("halfmax")
halfmax$halfmax= halfmax$halfmax %>% as.character %>% as.numeric

#cast back into wide format based on isotype
wide.halfmax=dcast(halfmax, RIX + dam + sire + ID + day + virus + antigen_virus + antigen + cohort + assay_date ~ isotype,mean,value.var='halfmax')
write.csv(wide.halfmax,"data_bin/halfmax_by_isotype.csv",row.names=F)
write.csv(halfmax,"data_bin/halfmax_data.csv",row.names=F)

# #### look for messed up duplicates ####
# dup=which(duplicated(lowestdil[c("RIX","ID","day","isotype")]) | 
#             duplicated(lowestdil[c("RIX","ID","day","isotype")], fromLast = TRUE))
# lowestdil.dup=lowestdil[dup,]
# lowestdil.dup=lowestdil.dup[order(lowestdil.dup$RIX,lowestdil.dup$day,lowestdil.dup$isotype,lowestdil.dup$ID),]
# write.csv(lowestdil.dup,"halfmax_duplicated.csv",row.names=F)


#### LAST POSITIVE (absorbance > 0.2) ####

#function takes anything under the threshold and sets it to 5 (arbitrary high number)
#then uses which.min to pick out the index of the lowest value in that column
#note: outputs INDEX of row, not actual log dilution value
dilfind=function(x){x[x<0.2] = 5;which.min(x)}

#apply over each column of the data frame
dil=NULL
dil=apply(data1[(length(data1)-7):length(data1)],1,dilfind)
dil=as.numeric(dil)

#if the row with the lowest dilution has a value less than threshold, replace that row index with 0
#using "0" instead of NA to distinguish between sample not run
for (i in 1:nrow(data1)){
  if(data1[i,(dil[[i]]+hc)]<0.2)
  {dil[[i]]=0}}

#create new lastpos data frame
lastpos=cbind(data1[1:hc],dil)
colnames(lastpos)[(ncol(lastpos))]="lastpos"

#convert lastpos value from column index to dilution factor
colnames(dilframe)[1]=c("lastpos")
lastpos=merge(lastpos,dilframe,all=T)
lastpos=lastpos[-1]
colnames(lastpos)[length(lastpos)]=c("lastpos")
lastpos$lastpos= lastpos$lastpos %>% as.character %>% as.numeric


#cast back into wide format based on isotype
wide.lastpos=dcast(lastpos, RIX + dam + sire + ID + day + virus + antigen_virus + antigen + cohort + assay_date ~ isotype,mean,value.var='lastpos')
write.csv(wide.lastpos,"data_bin/last_positive_by_isotype.csv",row.names=F)
write.csv(lastpos,"data_bin/lastpos_data.csv",row.names=F)

# #### look for messed up duplicates ####
# dup=which(duplicated(lowestdil[c("RIX","ID","day","isotype")]) | 
#             duplicated(lowestdil[c("RIX","ID","day","isotype")], fromLast = TRUE))
# lowestdil.dup=lowestdil[dup,]
# lowestdil.dup=lowestdil.dup[order(lowestdil.dup$RIX,lowestdil.dup$day,lowestdil.dup$isotype,lowestdil.dup$ID),]
# write.csv(lowestdil.dup,"lastpos_duplicated.csv",row.names=F)

# maybe useful to have merged last pos and lp data together?
# hm.lp=merge(halfmax,lastpos)

### cleanup ###
# rm(CC_names,data1,data1.nomocks,infections,RIX_CC,
#    RIX_CC_1,RIX_CC_2,RIX_names,data1_CC,newrow,RIX_CC_names)
# # rm(aucdata.dup,lowestdil.dup,halfmax.dup,dup,dup.auc)
# rm(dilutions,data1)
