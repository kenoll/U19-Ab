library(dplyr)
library(doBy)

#call up master founder mosaic file
setwd("~/Dropbox/Heise/U19-Ab/qtl_info/") 
mosaics=read.csv("~/Dropbox/Heise/CC/cc_founder_mosaics/founder_mosaics.csv")
mosaics$strain=gsub("/.+","",mosaics$strain)

#functions
region.haplos=function(chromo,start,end){
  region.2=subset(mosaics,
                  mosaics$start_position<=start & 
                    mosaics$end_position>=start & 
                    mosaics$chromosome==chromo)
  
  region.3=subset(mosaics,
                  mosaics$start_position<=end & 
                    mosaics$end_position>=end & 
                    mosaics$chromosome==chromo)
  
  region.4=subset(mosaics,
                  mosaics$start_position>=start & 
                    mosaics$end_position<=end & 
                    mosaics$chromosome==chromo)
  
  region=merge(region.2,region.3,all=T)
  region=merge(region,region.4,all=T)  
  return(region)
}


# 
get.haploscore=function(chromo,start,end,allele.effects){
  region=region.haplos(chromo,start,end)
  founders=c("A","B","C","D","E","F","G","H")
  CC.effects=data.frame(founders,allele.effects)
  
  #eliminate second haplotype if strain has same founder for both (most cases)
  #if there's a double recombination so that the same founders repeat, it will only keep one of each
  haplos=merge(region,CC.effects,by="founders")
  haplos=haplos[order(haplos$strain,haplos$chromosome,haplos$haplotype,haplos$start_position),]
  haplos=haplos[c(2:length(haplos),1)]
  haplos$strain_founder=paste(haplos$strain,haplos$founder,sep="_")
  haplos = haplos[!duplicated(haplos$strain_founder),]
  haplos = haplos[1:(length(haplos)-1)]
  return(haplos)
}

#if you only want to look at un-ambiguous calls (no recombinations, no heterozygosity)
get.haploscore2=function(chromo,start,end,allele.effects){
  haplos2=get.haploscore(chromo,start,end,allele.effects)
  haplos2=haplos2[!(duplicated(haplos2$strain) | duplicated(haplos2$strain, fromLast = TRUE)),]
  return(haplos2)
}

get.multi=function(chromo,start,end,allele.effects){
  haplos=get.haploscore(chromo,start,end,allele.effects)
  haplos2=get.haploscore2(chromo,start,end,allele.effects)
  multi=anti_join(haplos,haplos2)
  multi=multi[order(multi$strain,multi$chromosome,multi$haplotype,multi$start_position),]
  return(multi)
}

#### for effects only, not haplotypes ####
get.quickeffects=function(chromo,start,end,allele.effects){
  haplos=get.haploscore(chromo,start,end,allele.effects)
  score.avg=summaryBy(status~strain+alias+chromosome+haplotype,
                        data=haplos, FUN=mean, na.rm=T)
  score.avg=summaryBy(status.mean~strain+alias,
                        data=score.avg, FUN=mean, na.rm=T)
  colnames(score.avg)[length(score.avg)] = "status"
  return(score.avg)
}

