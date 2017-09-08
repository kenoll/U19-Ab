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
    trans.phenotype<-log(phenotype,base=base)
    dat.trans<-cbind(dat.trans, trans.phenotype)
  }
  
  return(dat.trans)
}
  

# hist.ab(df)
# for(i in 1:6)
# {
#   phenotype<-dat.inv[,i+hc]
#   hist(phenotype, main=abs[i])
# }  
#   
#   
#   colnames(dat.trans)<-col.list