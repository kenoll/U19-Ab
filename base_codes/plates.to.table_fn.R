plates.to.table=function(df){
  df = data.frame(df[,-1], row.names=df[,1])
  df[9,]=colnames(df)
  df=df[c(9,1:8),]
  df=as.data.frame(t(df))
  colnames(df)[1]="ID"
  return(df)
}