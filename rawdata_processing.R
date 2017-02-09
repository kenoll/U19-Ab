#### loads ####
library(dplyr)
library(reshape2)
library(tidyr)
library(flux)

#### getting started ####
#load file and make first column (dilutions) the row name
#the loaded data file has the first colum as the dilution vales
#the bottom row is the dates on which the ELISAs were run
#samples with wonky curves have been removed manually already (as per QC graphs at the end of the file)
setwd("~/Dropbox/Heise/ELISA Antibody/")
data1<-read.csv("cumulative_raw_ELISA_data_horizontal_bgsubtracted_curated_dates.csv")

#make dilutions row names
data1 = data.frame(data1[,-1], row.names=data1[,1])

#convert absorbance data to numeric values
data1[1:length(data1)] = sapply(data1[1:length(data1)], as.character)

#### Add infection status (Flu vs. Mock) and change RIX names to CC numbers ####
data1[10,]=colnames(data1)
data1=data1[c(10,9,1:8),]
data1=as.data.frame(t(data1))
colnames(data1)[1]="ID"

data1=separate(data1,ID,into=c("RIX","ID","day","isotype"),sep="_",extra="drop")

### infection status ###
#add virus infection status and remove mocks
#NOTE: deletes all rows for which there is no data about infection status

# infections=read.csv("~/Dropbox/Heise/weight_loss/weights_2016_10.csv")
# infections=infections[c(1,2,7,6,5,3)]
# colnames(infections)=c("RIX","ID","day","virus","cohort","RIX_ID")
# write.csv(infections,"infection_status.csv",row.names=F)
# #manually concatenate to add X in front of RIX to merge w/ ELISA data struction
infections=read.csv("infection_status.csv")

# infections=read.csv("infection_status.csv")
data1.inf=merge(data1,infections,all.x=T)
data1.inf=data1.inf[c(1:5,14,6:13)] #reorder
data1.inf=subset(data1.inf,data1.inf$RIX!="AntiCA04" & data1.inf$RIX!="B6xB6")

##pull out rows w/ no attached infection status
data1.novirus=data1.inf[is.na(data1.inf$virus),]
data1.novirus=subset(data1.novirus,data1.novirus$isotype=="IgG2ac")
write.csv(data1.novirus, "no_virus.csv",row.names=F)

data1.nomocks=subset(data1.inf,data1.inf$virus=="Influenza")
data1=data1.nomocks

#### import CC names ####
CC_names=read.csv("~/Dropbox/Heise/CC/cc_names.csv")
CC_names$Alias=as.character(CC_names$Alias)
CC_names$CCLine=as.character(CC_names$CCLine)

RIX_names <- data.frame(do.call("rbind", strsplit(data1$RIX,"x")))
RIX_names[,1] = gsub("X","",RIX_names[,1])
RIX_names[,1]=as.character(RIX_names[,1])
RIX_names[,2]=as.character(RIX_names[,2])

RIX_CC_1=data.frame(RI_1=RIX_names$X1, CC_1=CC_names[match(RIX_names$X1, CC_names$Alias),1])
RIX_CC_2=data.frame(RI_2=RIX_names$X2, CC_2=CC_names[match(RIX_names$X2, CC_names$Alias),1])

RIX_CC=cbind(RIX_CC_1,RIX_CC_2)
RIX_CC_names=paste(RIX_CC[,2],RIX_CC[,4],sep="x")

data1[,1]=RIX_CC_names
data1_CC=as.vector(paste(data1[,1],data1[,2],data1[,3],data1[,4],sep="_"))

#make sure all are numeric
data1[7:14] = sapply(data1[7:14], as.character)
data1[7:14] = sapply(data1[7:14], as.numeric)

#### Calculate AUC ####

#make vector for dilutions to do calculations with later
dilutions=c(1/100,1/300,1/1000,1/3000,1/10000,1/30000,1/100000,1/300000)
dilutions=abs(log10(dilutions))

#make a new row in the df to populate with AUC data
newrow = 0
data1 = cbind(data1,newrow)

#calculate AUC and add to data frame in last row
for(i in 1:nrow(data1))
{
  data1[i,15]<-auc(dilutions, data1[i,7:14])
}  

#create new DF with sample IDs and AUC data only and export
colnames(data1)[15]="AUC"
aucdata=data1[c(1:6,15)]
write.csv(aucdata,"AUC_data.csv",row.names=F)

#remove AUC from data
data1=data1[1:14]

#pull out duplicates to examine for QC
dup.auc=which(duplicated(aucdata[c("RIX","ID","day","isotype")]) | 
                duplicated(aucdata[c("RIX","ID","day","isotype")], fromLast = TRUE))
aucdata.dup=aucdata[dup.auc,]
aucdata.dup=aucdata.dup[order(aucdata.dup$RIX,aucdata.dup$ID,aucdata.dup$day,aucdata.dup$isotype),]

write.csv(aucdata.dup,"AUC_duplicated.csv",row.names=F)

#cast back into wide format based on isotype
#will average out any duplicate values (e.g. same sample run different days)
wideauc=dcast(aucdata, RIX + ID + day + assay_date ~ isotype,mean,value.var='AUC')
write.csv(wideauc,"AUC_by_isotype.csv",row.names=F)

#### Last Positive Dilution (absorbance >1.5) ####
#function takes anything under the threshold and sets it to 5 (arbitrary high number)
#then uses which.min to pick out the index of the lowest value in that column
#note: outputs INDEX of row, not actual log dilution value
dilfind=function(x){x[x<1.75] = 5;which.min(x)}

#apply over each column of the data frame
lowestdil=NULL
lowestdil=apply(data1[7:14],1,dilfind)

#if the row with the lowest dilution has a value less than threshold, replace that row index with 0
#using "0" instead of NA to distinguish between sample not run
for (i in 1:nrow(data1)){
  if(data1[i,(lowestdil[[i]]+6)]<1.75)
  {lowestdil[[i]]=0}}

#export data
data1[15]=lowestdil
colnames(data1)[15]="halfmax"
halfmax=data1[c(1:6,15)]
data1=data1[-15]

#cast back into wide format based on isotype
wide.halfmax=dcast(halfmax, RIX + ID + day + assay_date + virus ~ isotype,mean,value.var='halfmax')
wide.halfmax=wide.halfmax[c(1:4,6:11)]
write.csv(wide.halfmax,"halfmax_by_isotype.csv",row.names=F)

# #### look for messed up duplicates ####
# dup=which(duplicated(lowestdil[c("RIX","ID","day","isotype")]) | 
#             duplicated(lowestdil[c("RIX","ID","day","isotype")], fromLast = TRUE))
# lowestdil.dup=lowestdil[dup,]
# lowestdil.dup=lowestdil.dup[order(lowestdil.dup$RIX,lowestdil.dup$day,lowestdil.dup$isotype,lowestdil.dup$ID),]
# write.csv(lowestdil.dup,"halfmax_duplicated.csv",row.names=F)



#### LAST POSITIVE ####

#function takes anything under the threshold and sets it to 5 (arbitrary high number)
#then uses which.min to pick out the index of the lowest value in that column
#note: outputs INDEX of row, not actual log dilution value
dilfind=function(x){x[x<0.2] = 5;which.min(x)}

#apply over each column of the data frame
lowestdil=NULL
lowestdil=apply(data1[7:14],1,dilfind)

#if the row with the lowest dilution has a value less than threshold, replace that row index with 0
#using "0" instead of NA to distinguish between sample not run
for (i in 1:nrow(data1)){
  if(data1[i,(lowestdil[[i]]+6)]<0.2)
  {lowestdil[[i]]=0}}

#export data
data1[15]=lowestdil
colnames(data1)[15]="last_positive"
lastpos=data1[c(1:6,15)]
data1=data1[-15]

#cast back into wide format based on isotype
wide.lastpos=dcast(lastpos, RIX + ID + day + assay_date + virus ~ isotype,mean,value.var='last_positive')
wide.lastpos=wide.lastpos[c(1:4,6:11)]
write.csv(wide.lastpos,"last_positive_by_isotype.csv",row.names=F)

# #### look for messed up duplicates ####
# dup=which(duplicated(lowestdil[c("RIX","ID","day","isotype")]) | 
#             duplicated(lowestdil[c("RIX","ID","day","isotype")], fromLast = TRUE))
# lowestdil.dup=lowestdil[dup,]
# lowestdil.dup=lowestdil.dup[order(lowestdil.dup$RIX,lowestdil.dup$day,lowestdil.dup$isotype,lowestdil.dup$ID),]
# write.csv(lowestdil.dup,"lastpos_duplicated.csv",row.names=F)


### cleanup ###
rm(aucdata.dup,CC_names,data1.inf,data1.nomocks,infections,lowestdil.dup,RIX_CC,
   RIX_CC_1,RIX_CC_2,RIX_names,data1_CC,dup,dup.auc,newrow,RIX_CC_names)
rm(dilutions,data1,data1.novirus)

#### save important objects ####
save(wideauc,file="wideauc.Rdata")
save(wide.lastpos,file="wide.lastpos.Rdata")
save(wide.halfmax,file="wide.halfmax.Rdata")