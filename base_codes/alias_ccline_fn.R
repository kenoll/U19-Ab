CC_names=read.csv("~/Dropbox/Heise/CC/cc_names_a2l.csv")
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
  
  df["RIX"]=RIX_CC$RIX
  
  return(df)
}

alias.to.line.base=function(df){
  df$dam=as.character(df$dam)
  df$sire=as.character(df$sire)
  
  colnames(CC_names)[1]="sire"
  df=merge(df,CC_names[1:2],by="sire",all.x=T)
  df["sire"]=df["CCLine"]
  df=df[-length(df)]
  
  colnames(CC_names)[1]="dam"
  df=merge(df,CC_names[1:2],by="dam",all.x=T)
  df["dam"]=df["CCLine"]
  df=df[-length(df)]
  
  df["RIX"]=paste(df$dam,df$sire,sep="x")
  return(df)
}


alias.to.line2=function(df){
  if(any("RIX" %in% colnames(df)) && any(!("dam" %in% colnames(df)))) {
    df=split.rix(df)
    df=alias.to.line.base(df)
    df=df[-(1:2)]
  } else {
    if(any("dam" %in% colnames(df))){
      df=alias.to.line.base(df)
    }
  }
  return(df)
}
  


CC_names2=read.csv("~/Dropbox/Heise/CC/cc_names.csv")
CC_names2[1:2] = sapply(CC_names2[1:2], as.character)

line.to.alias=function(df){
  df$dam=as.character(df$dam)
  df$sire=as.character(df$sire)
  
  colnames(CC_names2)[1]="sire"
  df=merge(df,CC_names2[1:2],by="sire",all.x=T)
  df["sire"]=df["Alias"]
  df=df[-length(df)]
  
  colnames(CC_names2)[1]="dam"
  df=merge(df,CC_names2[1:2],by="dam",all.x=T)
  df["dam"]=df["Alias"]
  df=df[-length(df)]
  
  df["RIX"]=paste(df$dam,df$sire,sep="x")
  
  return(df)
}