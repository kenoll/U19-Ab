ccstatus=read.csv("~/Dropbox/Heise/CC/CCstatus.csv")

colnames(ccstatus)[6]="strain"
ccstatus$strain=sub("/.*","",ccstatus$strain)

ccstatus=ccstatus[c(2:3,6,27)]


setwd("/Users/kelseynoll/Dropbox/Heise/U19-Ab/weight/qtl_info/")
status=read.csv("chr1_mx1_combo_scores.csv")

status=merge(status,ccstatus,all.x=T)
write.csv(status,"chr1_recomb_mx1_combo_scores_avail.csv",row.names=F)
