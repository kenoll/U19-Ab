CC_names=read.csv("~/Dropbox/Heise/CC/cc_names.csv")
  CC_names[1:2] = sapply(CC_names[1:2], as.character)

alias.to.line=function(df){
  df$RIX=as.character(df$RIX)
  RIX_names <- data.frame(do.call("rbind", strsplit(df$RIX,"x")))
  RIX_names[,1] = gsub("X","",RIX_names[,1])
  RIX_names[1:2] = sapply(RIX_names[1:2], as.character)
  
  RIX_CC_1=data.frame(RI_1=RIX_names$X1, dam=CC_names[match(RIX_names$X1, CC_names$Alias),1])
  RIX_CC_2=data.frame(RI_2=RIX_names$X2, sire=CC_names[match(RIX_names$X2, CC_names$Alias),1])
  
  RIX_CC=cbind(RIX_CC_1,RIX_CC_2)
  RIX_CC$RIX=paste(RIX_CC[,2],RIX_CC[,4],sep="x")
  
  df[,1]=RIX_CC$RIX
  
  return(df)
}
