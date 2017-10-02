plates.to.table=function(df,date=F){
  df = data.frame(df[,-1], row.names=df[,1])
 
  if (date==F){
    df[9,]=colnames(df)
    df=df[c(9,1:8),]
    df=as.data.frame(t(df))
    colnames(df)[1]="ID"
    
    df[2:length(df)] = sapply(df[2:length(df)], as.character)
    df[2:length(df)] = sapply(df[2:length(df)], as.numeric)
  }
  
  if (date==T){
    df[10,]=colnames(df)
    df=df[c(10,9,1:8),]
    df=as.data.frame(t(df))
    colnames(df)[1]="ID"
    
    df[3:length(df)] = sapply(df[3:length(df)], as.character)
    df[3:length(df)] = sapply(df[3:length(df)], as.numeric)
  } 
  return(df)
}



