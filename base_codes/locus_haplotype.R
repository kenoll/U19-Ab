library(dplyr)
library(doBy)

#call up master founder mosaic file
setwd("~/Dropbox/Heise/CC/cc_founder_mosaics")
mosaics=read.csv("founder_mosaics.csv")

mosaics$strain=gsub("/.+","",mosaics$strain)

#selection region of interest based on chromosome number, start position, and end position
#will not include strain if the region is not in a 'founder block' i.e. spans multiple rows of the data frame
# chrom=1
# start=38509260
# end=40243410

chrom=15
start=75044018
end=75048853

directory="~/Dropbox/Heise/U19-Ab/weight/qtl_info/"

# set regions
region.1=subset(mosaics,
                mosaics$start_position<=start &
                  mosaics$end_position>=end &
                  mosaics$chromosome==chrom)

#capture founder blocks that overlap with region of interest anywhere
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

rm(region.2,region.3,region.4)


### select status of region based on founder by uncommenting applicable code (two or three ways)

# for two cases only (e.g. positive for staining), change/add/delete letters for founder strains to be marked as 'true'
group1.list=c("D","F","G")
group1="0"

group2.list=c("A","B","C","E","H")
group2="1"

region$status=
  ifelse(region$founder %in% group1.list
         ,group1,group2)

 
# ## for a three way case (e.g. mx1 status)
# group1="1"
# group2="0"
# group3="2"
# 
# group1.list=c("B","D","G")
# group1="2"
# 
# group2.list=c("A","C","H")
# group2="1"
# 
# group3.list=c("F","E")
# group3="0"
# 
# region$status=
#   ifelse(region$founder %in% group1.list
#          ,group1,
#          ifelse(region$founder %in% group2.list
#                 ,group2,group3)
#   )


#### to give founder haplotypes and allele score at locus ####

#eliminate second haplotype if strain has same founder for both (most cases)
#if there's a double recombination so that the same founders repeat, it will only keep one of each
haplos=region
haplos$alias_founder=paste(haplos$alias,haplos$founder,sep="_")
haplos = haplos[!duplicated(haplos$alias_founder),]
haplos = haplos[1:(length(haplos)-1)]

write.csv(haplos,file=paste0(directory,"chr",chrom,"_",
                               floor(start/1000000),"-",ceiling(end/1000000),"_haplotypes.csv"),row.names=F)

#if you only want to look at un-ambiguous calls (no recombinations, no heterozygosity)
haplos.2=haplos
haplos.2=haplos.2[!(duplicated(haplos.2$strain) | duplicated(haplos.2$strain, fromLast = TRUE)), ]
haplos.2=haplos.2[c(1,2,7,8)]

write.csv(haplos.2,file=paste0(directory,"chr",chrom,"_",
                                 floor(start/1000000),"-",ceiling(end/1000000),"_norecomb_haplotypes.csv"),row.names=F)

`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

eliminated.strains=subset(haplos$strain,haplos$strain %not in% haplos.2$strain)

library(dplyr) 
dropped=anti_join(haplos,haplos.2)
write.csv(dropped,
          file=paste0(directory,"chr",chrom,"_",
                      floor(start/1000000),"-",ceiling(end/1000000),
                      "_recomb_haplotypes.csv"),row.names=F)


## for specific strains of interest only
# strain.list=c("CC030","CC044","CC031","CC011","CCO75")
# 
# haplos.strains=haplos[which(haplos$strain %in% strain.list),]
# haplos.strains
# write.csv(haplos.strains,"~/Desktop/haplos.csv")
# region.strains=region[which(region$strain %in% strain.list),]
# region.strains

#### for scores only, not haplotypes ####
library(doBy)
region$status=as.numeric(region$status)

region.avgs=summaryBy(status~strain+alias+chromosome+haplotype,
                      data=region, FUN=mean, na.rm=T)
region.avgs=summaryBy(status.mean~strain+alias,
                      data=region.avgs, FUN=mean, na.rm=T)
colnames(region.avgs)[length(region.avgs)] = "status"

region.avgs$strain=as.character(region.avgs$strain)

head(region.avgs)

directory="~/Dropbox/Heise/weight_loss/qtls/"
write.csv(region.avgs,file=paste(directory,"chr",chrom,"_",
          floor(start/1000000),"-",ceiling(end/1000000),"_scores.csv",sep=""),row.names=F)
