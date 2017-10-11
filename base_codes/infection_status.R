setwd("~/Dropbox/Heise/U19-Ab")
infections1=read.csv("weight/weights_2017_05.csv")
infections1=infections1[c(1,2,7,6,5,3)]
infections2=alias.to.line2(infections1)

colnames(infections1)=c("Alias.RIX","ID","day","virus","cohort","RIX_ID")
colnames(infections2)=c("RIX","ID","day","virus","cohort","RIX_ID")

infections=merge(infections2,infections1)
infections=infections[order(infections$RIX,infections$day,infections$ID),]
infections=infections[c(6,7,5,1:4)]
# infections$RIX=paste0("X",infections$RIX)
write.csv(infections,"antibody/data_bin/infection_status.csv",row.names=F)