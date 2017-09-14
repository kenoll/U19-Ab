is.zero=function(x){
  ifelse(x==0,T,F)
}


trans.ab=function(df,exp,hc=hc,num=6){
  dat.trans=df[,1:hc]

  for (i in 1:num)
  {
    phenotype<-df[,i+hc]
    trans.phenotype<-(phenotype)^exp
    dat.trans<-cbind(dat.trans, trans.phenotype)
  }
  
  return(dat.trans)
}


###

log.ab=function(df,hc=hc,base=exp(1),num=6){
  dat.trans=df[,1:hc]
  
  for (i in 1:num)
  {
    phenotype<-df[,i+hc]
    phenotype[phenotype<=0] <- min(subset(phenotype,phenotype>0),na.rm=T)
    trans.phenotype<-log(phenotype,base=base)
    dat.trans<-cbind(dat.trans, trans.phenotype)
  }
  
  return(dat.trans)
}
  
###
hist.ab=function(df,hc=hc,num=6){
  for(i in 1:num)
  {
    phenotype<-df[,i+hc]
    hist(phenotype, main=paste(gsub("dat*","",deparse(substitute(df))),abs[i],sep=" "))
  }
}



# ###
# box.ab=function(df,hc=hc,num=6){
#   for (i in 1:num){
#     phenotype<-df[,i+hc]
#     phenotype[phenotype<=0] <- min(subset(phenotype,phenotype>0),na.rm=T)
#     lm.bc=lm(phenotype~RIX,data=df)
#     bc=boxcox(lm.bc)
#     print(bc$x[which.max(bc$y)])
#     
#   }
# }