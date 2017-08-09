mx1=read.csv("~/Dropbox/Heise/CC/mx1_status.csv")

split.rix=function(df) {
  if(any(!("RIX" %in% colnames(df)))) stop('RIX is not a column name in the dataset')
  df$RIX=as.character(df$RIX)
  RIX_sep <- data.frame(do.call("rbind", strsplit(df$RIX,"x")))
  colnames(RIX_sep)[1:2]=c("dam","sire")
  df=cbind(RIX_sep,df)
  return(df)
}

add.stat.base=function(df,locus.score,locus.name) {
  if(any("founders" %in% colnames(locus.score))){
    locus.score=locus.score[c(colnames(locus.score[1]),"status","founders")]
    
    colnames(locus.score)[1]="sire"
    df=merge(df,locus.score)
    names(df)[names(df) == 'status']=paste0("sire_",locus.name)
    names(df)[names(df) == 'founders']=paste0("sire.founder_",locus.name)
    
    colnames(locus.score)[1]="dam"
    df=merge(df,locus.score)
    names(df)[names(df) == 'status']=paste0("dam_",locus.name)
    names(df)[names(df) == 'founders']=paste0("dam.founder_",locus.name)
  } else {
    locus.score=locus.score[c(colnames(locus.score[1]),"status")]
    
    colnames(locus.score)[1]="sire"
    df=merge(df,locus.score)
    names(df)[names(df) == 'status']=paste0("sire_",locus.name)
    
    colnames(locus.score)[1]="dam"
    df=merge(df,locus.score)
    names(df)[names(df) == 'status']=paste0("dam_",locus.name)
  }
  
  df$additive=rowSums(df[c(paste0("dam_",locus.name),paste0("sire_",locus.name))])
  df$additive=as.character(df$additive)
  df$additive=as.factor(df$additive)
  names(df)[names(df) == 'additive']=paste0("add_",locus.name)
  
  labels <- apply(df[,c(paste0("dam_",locus.name),paste0("sire_",locus.name))], 1, sort)
  df$combo<- factor(apply(labels, 2, function(x) paste(x, collapse="_")))
  names(df)[names(df) == 'combo']=paste0("combo_",locus.name)
  
  return(df)
  }
  

add.stat.gene=function(df,locus.score,locus.name) {
  if(any("RIX" %in% colnames(df)) && any(!("dam" %in% colnames(df)))) {
    df=split.rix(df)
    df=add.stat.base(df,locus.score,locus.name)
  } else {
    if(any("dam" %in% colnames(df))){
      df=add.stat.base(df,locus.score,locus.name)
    } else {
      colnames(locus.score)[1]=colnames(df)[1]
      df=merge(df,locus.score)
      names(df)[names(df) == 'status']=locus.name
    }
  }
  return(df)
}


add.stat.mx1 = function(df,locus.score,locus.name,cast.score=0.5){
  if (cast.score==1) cast.col=9
  if (cast.score==0.5) cast.col=10
  if (cast.score==0) cast.col=11
  if (missing(cast.score)) cast.col=10
  
  mx1.score=mx1[c(1,2,cast.col)]
  colnames(mx1.score)[3]="status"

  add.stat.gene(df,locus.score=mx1.score,locus.name="mx1")
}

