library(dplyr)
library(doBy)

#call up master founder mosaic file
setwd("~/Dropbox/Heise/U19-Ab/qtl_info/") 
mosaics=read.csv("~/Dropbox/Heise/CC/cc_founder_mosaics/founder_mosaics.csv")
mosaics$strain=gsub("/.+","",mosaics$strain)

# set regions of interest
chrom=1
start=38509260
end=40243410
locus.name="chr1"

#strain scores
strain.scores=c(
A.score=1 ,
B.score=2 ,
C.score=1 ,
D.score=2 ,
E.score=0 ,
F.score=0 ,
G.score=2 ,
H.score=1
)

get.regions=function(chrom,start,end){
  region.2=subset(mosaics,
                  mosaics$start_position<=start & 
                    mosaics$end_position>=start & 
                    mosaics$chromosome==chrom)
  
  region.3=subset(mosaics,
                  mosaics$start_position<=end & 
                    mosaics$end_position>=end & 
                    mosaics$chromosome==chrom)
  
  region.4=subset(mosaics,
                  mosaics$start_position>=start & 
                    mosaics$end_position<=end & 
                    mosaics$chromosome==chrom)
  
  region=merge(region.2,region.3,all=T)
  region=merge(region,region.4,all=T)  
  return(region)
}

# 
get.haplos=function(chrom,start,end){
  region=get.regions(chrom,start,end)
  founders=c("A","B","C","D","E","F","G","H")
  status=strain.scores
  CC.scores=data.frame(founders,status)
  
  #eliminate second haplotype if strain has same founder for both (most cases)
  #if there's a double recombination so that the same founders repeat, it will only keep one of each
  haplos=merge(region,CC.scores,by="founders")
  haplos=haplos[order(haplos$strain,haplos$chromosome,haplos$haplotype,haplos$start_position),]
  haplos=haplos[c(2:length(haplos),1)]
  haplos$strain_founder=paste(haplos$strain,haplos$founder,sep="_")
  haplos = haplos[!duplicated(haplos$strain_founder),]
  haplos = haplos[1:(length(haplos)-1)]
  return(haplos)
}

#if you only want to look at un-ambiguous calls (no recombinations, no heterozygosity)
get.haplos2=function(chrom,start,end){
  haplos2=get.haplos(chrom,start,end)
  haplos2=haplos2[!(duplicated(haplos2$strain) | duplicated(haplos2$strain, fromLast = TRUE)),]
  return(haplos2)
}

get.multi=function(chrom,start,end){
  haplos=get.haplos(chrom,start,end)
  haplos2=get.haplos2(chrom,start,end)
  multi=anti_join(haplos,haplos2)
  multi=multi[order(multi$strain,multi$chromosome,multi$haplotype,multi$start_position),]
  return(multi)
}

#### for scores only, not haplotypes ####
get.quickscores=function(chrom,start,end){
  haplos=get.haplos(chrom,start,end)
  score.avg=summaryBy(status~strain+alias+chromosome+haplotype,
                        data=haplos, FUN=mean, na.rm=T)
  score.avg=summaryBy(status.mean~strain+alias,
                        data=score.avg, FUN=mean, na.rm=T)
  colnames(score.avg)[length(score.avg)] = "status"
  return(score.avg)
}