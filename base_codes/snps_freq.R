library(ggplot2)
library(plyr)
# library(dplyr)

#sanger downloaded file w/ strain IDs added to Csq columns - either all mutations or missense only

snps=read.csv("~/Dropbox/Heise/ELISA Antibody/qtls/d10/pre-feb2017/snps/snps_11-69082450-72556480.csv"
              , na.strings=c("","-","~"," "))

#OR

snps=read.csv("~/Dropbox/Heise/ELISA Antibody/qtls/d10/pre-feb2017/snps/snps_11-69082450-72556480_missenseonly.csv"
              , na.strings=c("","-","~"," "))

colnames(snps)=c("Chr"        ,  "Position"   ,  "Gene"       ,  "dbSNP"     ,   "Ref"         ,
                 "129S1_SvImJ",  "129_Csq"    ,  "A_J"        ,  "AJ_Csq"    ,   "CAST_EiJ"    ,
                 "CAST_Csq"   ,  "NOD_ShiLtJ" ,  "NOD_Csq"    ,  "NZO_HlLtJ" ,   "NZO_Csq"     ,
                 "PWK_PhJ"    ,  "PWK_Csq"    ,  "WSB_EiJ"    ,  "WSB_Csq"     )

#genes for which WSB has a private mutation
gene.list= c("Wscd1","Zbtb4","Zfp3","Tnk1","Rpain","Rpain","Slc16a13","Kctd11","Nlrp1b","Nlrp1a")

#frequency table (plyr count function works here, dlpyr count function does not)
snps.count=count(snps$Gene)
  colnames(snps.count)=c("Gene","Freq")

#compare frequency of mutations across region
barplot(snps.count$Freq, space=0, names = snps.count$Gene, xlab='Gene', ylab = 'Freq')

#pull out highly polymorphic genes
snps.count[which.max(snps.count$Freq),]
median(snps.count$Freq)
subset(snps.count,snps.count$Freq>5*median(snps.count$Freq))

#frequency of mutations across all strains only in subset of gene.list
snps.subset=snps[which(snps$Gene %in% gene.list),]
  snps.count=count(snps.subset$Gene)
    colnames(snps.count)=c("Gene","Freq")
  snps.count

#frequency of mutations in WSB (and any other strain)
snps.wsb=snps[complete.cases(snps["WSB_Csq"]),]
  # snps.wsb.count=count(snps.wsb$Gene)
  #   colnames(snps.wsb.count)=c("Gene","Freq")
  # snps.wsb.count
  
  #frequency of mutations in WSB (and any other strain) in gene.list
  snps.wsb.subset=snps.wsb[which(snps.wsb$Gene %in% gene.list),]
    snps.wsb.count=count(snps.wsb.subset$Gene)
      colnames(snps.wsb.count)=c("Gene","Freq")
    snps.wsb.count
  
  
#frequency of private WSB mutations
snps.wsb.private=snps.wsb[is.na(snps.wsb[6]) & is.na(snps.wsb[8]) & is.na(snps.wsb[10]) 
                          & is.na(snps.wsb[12]) & is.na(snps.wsb[14]) & is.na(snps.wsb[16]),]
  # snps.wsb.private.count=count(snps.wsb.private$Gene)
  #   colnames(snps.wsb.private.count)=c("Gene","Freq")
  # snps.wsb.private.count

  #frequency of private WSB mutations in gene.list
  snps.wsb.private.subset=snps.wsb.private[which(snps.wsb.private$Gene %in% gene.list),]
    snps.wsb.private.count=count(snps.wsb.private.subset$Gene)
      colnames(snps.wsb.private.count)=c("Gene","Freq")
    snps.wsb.private.counts
    
    
###############
    
    
#for multiple allele effect groups, testing w/ mx1 QTL region
snps=read.csv("~/Dropbox/Heise/weight_loss/sanger_snps_16-97018120-98096140.csv"
                  , na.strings=c("","-","~"," "))

    colnames(snps)=c("Chr","Position","Gene","dbSNP","Ref",
                     "129","129_Csq","AJ","AJ_Csq","CAST",
                     "CAST_Csq","NOD","NOD_Csq","NZO","NZO_Csq",
                     "PWK","PWK_Csq","WSB","WSB_Csq")

indels=read.csv("~/Dropbox/Heise/weight_loss/sanger_indels_16-97018120-98096140.csv"
              , na.strings=c("","-","~"," ")) 

    colnames(indels)=c("Chr","Position","Gene","dbSNP","Ref",
                     "129","129_Csq","AJ","AJ_Csq","CAST",
                     "CAST_Csq","NOD","NOD_Csq","NZO","NZO_Csq",
                     "PWK","PWK_Csq","WSB","WSB_Csq")
    
muts=rbind(snps,indels)

#set groups of interest, group those w/ B6 accordingly
group.1=c("NZO","PWK")
group.2=c("CAST")
group.B6=c("129","NOD","AJ","WSB")

muts.group1=muts[is.na(muts["NZO"]) | is.na(muts["PWK"]),]
# snps.wsb.private.count=count(snps.wsb.private$Gene)
#   colnames(snps.wsb.private.count)=c("Gene","Freq")
# snps.wsb.private.count