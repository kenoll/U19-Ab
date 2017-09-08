raw.weights=read.csv("weight/weights_2016_10.csv")
  raw.weights=subset(weights,weights$Treatment!="")
  raw.weights=raw.weights[,c(1:3,5:7,9:24,40)]
  raw.weights[7:23] = sapply(raw.weights[7:23], as.character)
  raw.weights[7:23] = sapply(raw.weights[7:23], as.numeric)

raw.weights$D1per=raw.weights$D1/raw.weights$D0*100
  raw.weights$D2per=raw.weights$D2/raw.weights$D0*100
  raw.weights$D3per=raw.weights$D3/raw.weights$D0*100
  raw.weights$D4per=raw.weights$D4/raw.weights$D0*100
  raw.weights$D5per=raw.weights$D5/raw.weights$D0*100
  raw.weights$D6per=raw.weights$D6/raw.weights$D0*100
  raw.weights$D7per=raw.weights$D7/raw.weights$D0*100
  raw.weights$D8per=raw.weights$D8/raw.weights$D0*100
  raw.weights$D9per=raw.weights$D9/raw.weights$D0*100
  raw.weights$D10per=raw.weights$D10/raw.weights$D0*100
  raw.weights$D11per=raw.weights$D11/raw.weights$D0*100
  raw.weights$D12per=raw.weights$D12/raw.weights$D0*100
  raw.weights$D13per=raw.weights$D13/raw.weights$D0*100
  raw.weights$D14per=raw.weights$D14/raw.weights$D0*100
  raw.weights$D15per=raw.weights$D15/raw.weights$D0*100
  raw.weights$D45per=raw.weights$D45/raw.weights$D0*100
  
raw.weights=raw.weights[c(1,2,3,5,6,7,24:30)]

raw.weights = alias.to.line(raw.weights)
raw.weights = split.rix(raw.weights)

strain.list.0=c("CC012","CC015","CC017","CC021","CC022","CC024","CC040","CC065")
strain.list.3=c("CC002","CC036","CC037","CC045","CC072")
strain.list=c(strain.list.0,strain.list.3)

weights=subset(raw.weights,raw.weights$dam %in% strain.list | raw.weights$sire %in% strain.list)
  weights=melt(weights,id=c("dam","sire","RIX","ID","RIX_ID","Treatment","Day"),value.name="weight")
  colnames(weights)[(length(weights)-1)]="timepoint"
  
  ggplot(data=weights,aes(x=timepoint,y=weight,color=Treatment,group=ID))+geom_point(na.rm=T)+
    facet_wrap("RIX",ncol=6)+labs(x="Day",y="Weight (% of starting)")+
    stat_smooth(method=loess,aes(group=Treatment,fill=Treatment))+
    # scale_y_continuous(breaks=seq(65,120,5))+
    ylim(60,130)+
    theme(axis.text.x=element_text(angle=60,vjust=1,hjust=1))+
    geom_hline(yintercept=100)+    
    ggtitle("Dam or Sire in Chr1 validation strain list")
  
  ggsave("qtl_info/chr1_validation_strains_weight.pdf",width=10,height=10)
  
sire.weights=subset(raw.weights,raw.weights$sire %in% strain.list)
  sire.weights=melt(sire.weights,id=c("dam","sire","RIX","ID","RIX_ID","Treatment","Day"),value.name="weight")
  colnames(sire.weights)[(length(sire.weights)-1)]="timepoint"
  
  ggplot(data=sire.weights,aes(x=timepoint,y=weight,color=Treatment,group=ID))+geom_point(na.rm=T)+
    facet_wrap("RIX",ncol=4)+labs(x="Day",y="Weight (% of starting)")+
    stat_smooth(method=loess,aes(group=Treatment,fill=Treatment))+
    # scale_y_continuous(breaks=seq(65,120,5))+
    ylim(60,130)+
    theme(axis.text.x=element_text(angle=60,vjust=1,hjust=1))+
    geom_hline(yintercept=100)+    
    ggtitle("Sires in chr1 validation strain list")


dam.weights=subset(raw.weights,raw.weights$dam %in% strain.list)
  dam.weights=melt(dam.weights,id=c("dam","sire","RIX","ID","RIX_ID","Treatment","Day"),value.name="weight")
  colnames(dam.weights)[(length(dam.weights)-1)]="timepoint"
  
  ggplot(data=dam.weights,aes(x=timepoint,y=weight,color=Treatment,group=ID))+geom_point(na.rm=T)+
    facet_wrap("RIX",ncol=4)+labs(x="Day",y="Weight (% of starting)")+
    stat_smooth(method=loess,aes(group=Treatment,fill=Treatment))+
    # scale_y_continuous(breaks=seq(65,120,5))+
    ylim(60,130)+
    theme(axis.text.x=element_text(angle=60,vjust=1,hjust=1))+
    geom_hline(yintercept=100)+
    ggtitle("Dams in chr1 validation strain list")
  

cross.weights=subset(raw.weights,raw.weights$dam %in% strain.list & raw.weights$sire %in% strain.list)
  cross.weights=melt(cross.weights,id=c("dam","sire","RIX","ID","RIX_ID","Treatment","Day"),value.name="weight")
  colnames(cross.weights)[(length(cross.weights)-1)]="timepoint"
  
  ggplot(data=cross.weights,aes(x=timepoint,y=weight,color=Treatment,group=ID))+geom_point(na.rm=T)+
    facet_wrap("RIX",ncol=2)+labs(x="Day",y="Weight (% of starting)")+
    stat_smooth(method=loess,aes(group=Treatment,fill=Treatment))+
    # scale_y_continuous(breaks=seq(65,120,5))+
    ylim(60,130)+
    theme(axis.text.x=element_text(angle=60,vjust=1,hjust=1))+
    geom_hline(yintercept=100)+
    ggtitle("Dam and Sire in chr1 validation strain list")


strain.list.2=c("CC060","CC037")
  
weights.spec=subset(raw.weights,raw.weights$dam %in% strain.list.2 | raw.weights$sire %in% strain.list.2)
  weights.spec=melt(weights.spec,id=c("dam","sire","RIX","ID","RIX_ID","Treatment","Day"),value.name="weight")
  colnames(weights.spec)[(length(weights.spec)-1)]="timepoint"
  
  ggplot(data=weights.spec,aes(x=timepoint,y=weight,color=Treatment,group=ID))+geom_point(na.rm=T)+
    facet_wrap("RIX")+labs(x="Day",y="Weight (% of starting)")+
    stat_smooth(method=loess,aes(group=Treatment,fill=Treatment))+
    # scale_y_continuous(breaks=seq(65,120,5))+
    ylim(60,130)+
    theme(axis.text.x=element_text(angle=60,vjust=1,hjust=1))+
    geom_hline(yintercept=100)+    
    ggtitle(paste0("Dam or Sire is ",strain.list.2))
  