#### to give haplotypes ####

##notes: certain user inputs are required and are described here
#line 7: set working directory/indicate where you have saved the master founder mosaic file
#line 11: chromosome number, start, and end location information for gene/locus/region of interest
#line 29: (optional) - see intructions for sorting out founders into categories (e.g. mx1 status or positive staining for a marker)
#line 55: indication of where you want to name/output the file to be generated


#call up master founder mosaic file
setwd("~/Dropbox/Heise/CC/cc_founder_mosaics")
mosaics=read.csv("founder_mosaics.csv")

#selection region of interest based on chromosome number, start position, and end position
#will not include strain if the region is not in a 'founder block' i.e. spans multiple rows of the data frame
chrom=11
# start=66376670
# end=66391310

start=70500000
end=72000000

region.1=subset(mosaics,
                mosaics$start_position<=start & 
                mosaics$end_position>=end & 
                mosaics$chromosome==chrom)

#to capture non-founder blocks
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
# region=merge(region.1,region.merge,all=T) #regions 2 and 3 already contain all in region 1

#region.diff highlights lines with a haplotype change within the region
library(dplyr)
# region.diff=setdiff(region,region.1)
region.diff=anti_join(region,region.1)
#will eliminate a non-changing second haplotype

#eliminate second haplotype if strain has same founder for both (most cases)
region$alias_founder=paste(region$alias,region$founder,sep="_")
# region = region[order(region$haplotype),] 
# region = region[order(region$strain),]
region = region[!duplicated(region$alias_founder),]
region = region[1:7]
#if there's a double recombination so that the same founders repeat, it will only keep one of each

region.diff$alias_founder=paste(region.diff$alias,region.diff$founder,sep="_")
region.diff = region.diff[!duplicated(region.diff$alias_founder),]
region.diff = region.diff[1:7]

## optional - select status of region based on founder by uncommenting applicable code (two or three ways)

# for two cases only (e.g. positive for staining), change/add/delete letters for founder strains to be marked as 'true'
region$status=
  ifelse(region$founder=="H",
         "0","1")

region.diff$status=
  ifelse(region.diff$founder=="H",
         "0","1")

# delete duplicates for lines where the founder haplotype change doesn't change the allele status
region.diff$alias_status=paste(region.diff$alias,region.diff$status,sep="_")
# region.diff = region.diff[!duplicated(region.diff$alias_status),]
# region.diff = region.diff[1:8]

changes=which(duplicated(region.diff["alias_status"])|(duplicated(region.diff["alias_status"],fromLast=T)))
changes=region.diff[changes,]
changes=setdiff(region.diff,changes)

# ## for a three way case (e.g. mx1 status)
# group1="dom"
# group2="mus"
# group3="cast"
# 
# region$status=
#   ifelse(region$founder=="A" |
#            region$founder=="B" |
#            region$founder=="C" |
#            region$founder=="D" |
#            region$founder=="H"
#           ,group1,
#       ifelse(region$founder=="E" |
#           region$founder=="G"
#           ,group2,group3)
#   )

#output csv
#user input where you want to save the file
file="~/Dropbox/Heise/ELISA Antibody/qtls/d10/chr5_62-75_status.csv"
write.csv(region,file,row.names=F)

file="~/Dropbox/Heise/ELISA Antibody/qtls/d10/chr5_transitionsonly_status.csv"
write.csv(region.diff,file,row.names=F)

file="~/Dropbox/Heise/ELISA Antibody/qtls/d10/chr5_meaningful_transitionsonly_status.csv"
write.csv(changes,file,row.names=F)







############## for scores only, not haplotypes ############

#call up master founder mosaic file
setwd("~/Dropbox/Heise/CC/cc_founder_mosaics")
mosaics=read.csv("founder_mosaics.csv")

#selection region of interest based on chromosome number, start position, and end position
#will not include strain if the region is not in a 'founder block' i.e. spans multiple rows of the data frame
chrom=11

start=70500000
end=72000000

region.1=subset(mosaics,
                mosaics$start_position<=start & 
                  mosaics$end_position>=end & 
                  mosaics$chromosome==chrom)

#to capture non-founder blocks
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
# region=merge(region.1,region.merge,all=T) #regions 2 and 3 already contain all in region 1

region$status=
  ifelse(region$founder=="H",
         "0","1")

region$status=as.numeric(region$status)

library(doBy)
region.avgs=summaryBy(status~strain+alias+chromosome+haplotype,
                      data=region, FUN=mean, na.rm=T)
region.avgs=summaryBy(status.mean~strain+alias,
                      data=region.avgs, FUN=mean, na.rm=T)
colnames(region.avgs)[length(region.avgs)] = "status"

region.avgs$strain=as.character(region.avgs$strain)


region.avgs$strain=gsub("/.+","",region.avgs$strain)
head(region.avgs)

write.csv(region.avgs,file=paste("~/Dropbox/Heise/ELISA Antibody/qtls/d10/snps/chr",chrom,"_",
          start/1000000,"-",end/1000000,"_scores.csv",sep=""),row.names=F)
