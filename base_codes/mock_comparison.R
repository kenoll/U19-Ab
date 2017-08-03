setwd("~/Dropbox/Heise/U19-Ab/antibody/phenotypes/mock/")

#read in allele score (as per founder haplotype) for locus of interest
chr.region="chr18_70-82"
chr.status=read.csv(file=paste("~/Dropbox/Heise/U19-Ab/antibody/data_bin/mock/",chr.region,
  "_scores.csv",sep=""))
chr.status=chr.status[c(1,3)]

#read in phenotype data
data.type="auc"
dat.7=read.csv(paste("~/Dropbox/Heise/U19-Ab/antibody/data_bin/d7/dat.7.",data.type,".csv",sep=""))
dat.10=read.csv(paste("~/Dropbox/Heise/U19-Ab/antibody/data_bin/d10/dat.10.",data.type,".csv",sep=""))
dat.15=read.csv(paste("~/Dropbox/Heise/U19-Ab/antibody/data_bin/d15/dat.15.",data.type,".csv",sep=""))
dat.45=read.csv(paste("~/Dropbox/Heise/U19-Ab/antibody/data_bin/d45/dat.45.",data.type,".csv",sep=""))

hc=7
abs=c("IgG1" ,  "IgG2ac" , "IgG2b" , "IgG3" ,  "IgM"  ,  "TotalG")
day.list=list(7,10,15,45)
dat.list=list(dat.7,dat.10,dat.15,dat.45)

b.e=NULL


#### Association Testing ####

for(k in 1:4)
{
  day=day.list[[k]]
  dat=dat.list[[k]]

  dat$RIX=as.character(dat$RIX)

  #associate genotype to phenotype
  dat$RIX=as.character(dat$RIX)
  RIX_sep <- data.frame(do.call("rbind", strsplit(dat$RIX,"x")))
  colnames(RIX_sep)[1:2]=c("dam","sire")
  
  dat=cbind(RIX_sep,dat)
  
  colnames(chr.status)[1]="dam"
  dat=merge(dat,chr.status)
  colnames(dat)[length(dat)]="dam.status"
  
  colnames(chr.status)[1]="sire"
  dat=merge(dat,chr.status)
  colnames(dat)[length(dat)]="sire.status"
  
  dat=dat[3:length(dat)]
  
  #inheritance models
  dat$additive.model=rowSums(dat[c("dam.status","sire.status")])
  
  for(i in 1:6)
  {
    lm.ab=lm(dat[,i+hc]~dat$additive.model)
    an=anova(lm.ab)
    print(an[1,5])
    b.eff=c(abs[i],an[1,5])
    b.e <- rbind(b.e, b.eff)
  }

}

  row.names(b.e)=1:nrow(b.e)
  b.e=as.data.frame(b.e)
  colnames(b.e)=c("isotype","p")
  b.e$day=rep(c("7","10","15","45"),each=length(abs))
  b.e=b.e[,c(3,1,2)]
  
  write.csv(b.e,file=(
    paste(chr.region,"_on_IAV_ab_",data.type,".csv",sep="")),
    row.names=F)


#### for setting significance threshold ####
#10x permutations
b.e=NULL
n.perms=100
  
for (j in 1:(n.perms+1))
{

for(k in 1:4)
{
  day=day.list[[k]]
  dat=dat.list[[k]]
  
  #associate genotype to phenotype
  dat$RIX=as.character(dat$RIX)
  RIX_sep <- data.frame(do.call("rbind", strsplit(dat$RIX,"x")))
  colnames(RIX_sep)[1:2]=c("dam","sire")
  
  dat=cbind(RIX_sep,dat)
  
  colnames(chr.status)[1]="dam"
  dat=merge(dat,chr.status)
  colnames(dat)[length(dat)]="dam.status"
  
  colnames(chr.status)[1]="sire"
  dat=merge(dat,chr.status)
  colnames(dat)[length(dat)]="sire.status"
  
  dat=dat[3:length(dat)]
  
  #inheritance models
  dat$additive.model=rowSums(dat[c("dam.status","sire.status")])
  
  #permute
  dat.perm=dat
  dat.perm$additive.model=sample(dat$additive.model)
  
  for(i in 1:6)
  {
    lm.ab=lm(dat[,i+hc]~dat.perm$additive.model)
    an=anova(lm.ab)
    # print(an[1,5])
    b.eff=c(abs[i],an[1,5])
    b.e <- rbind(b.e, b.eff)
  }
  
}
}
  
  row.names(b.e)=1:nrow(b.e)
  b.e=as.data.frame(b.e)
  colnames(b.e)=c("isotype","p")
  b.e$day=rep(c("7","10","15","45"),each=length(abs))
  b.e=b.e[,c(3,1,2)]
  b.e$iteration=rep(1:j,each=k*i)
  
  write.csv(b.e,file=(
  paste(chr.region,"_on_IAV_ab_",data.type,"_perms.csv",sep="")),
  row.names=F)


b.e=read.csv(file=(paste(chr.region,"_on_IAV_ab_",data.type,"_perms.csv",sep="")))
b.e.dat=read.csv(paste(chr.region,"_on_IAV_ab_",data.type,".csv",sep=""))
b.e.dat$iteration="NA"

quantile(b.e$p,probs=.05)

b.e=rbind(b.e,b.e.dat)

#plot distributions
bins=c(.0001,.0005,.001,.005,.01,.05,.1,.5,1)

# ggplot(b.e, aes(x=p)) + geom_histogram(bins=15,color="black",fill="gray75") + 
#   scale_x_log10(name="p-value",breaks=bins) + theme_bw() +
#   facet_wrap(~iteration,nrow=4) + theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1))

b.e.sub=b.e[1:1344,]

ggplot(b.e.sub, aes(x=p)) + geom_histogram(bins=15,color="black",fill="gray75") + 
  scale_x_log10(name="p-value",breaks=bins,limits=c(1e-5,1)) + theme_bw() +
  facet_wrap(~iteration,nrow=8) + theme(axis.text.x=element_text(angle=60,vjust=1,hjust=1))+
  ggtitle(paste0("p-value distributions for permutations of ",chr.region,"~ IAV-induced antibody across days/isotypes"))


#compare 

ggplot(b.e.dat, aes(x=p)) + geom_histogram(bins=15,color="black",fill="gray75") + 
  scale_x_log10(name="p-value",breaks=bins,limits=c(1e-5,1)) + theme_bw() +
  theme(axis.text.x=element_text(angle=60,vjust=1,hjust=1))+
  ggtitle(paste0("p-value distributions for ",chr.region,"~IAV-induced antibody across days/isotypes"))


##### mock data and pheno dist
    chr.region="chr18_70-82"
    chr.status=read.csv(file=paste("~/Dropbox/Heise/U19-Ab/antibody/data_bin/mock/",chr.region,
                               "_scores.csv",sep=""))
    colnames(chr.status)[1]="Strain"

    mock=read.csv("~/Dropbox/Heise/U19-Ab/antibody/data_bin/mock/mock_RI_dat.csv")
    colnames(mock)[1:3]=c("Strain","ID","Sex")
    colnames(mock)[4:12] = gsub("_ug","",colnames(mock[4:12]))
    colnames(mock)[5] = "TotalG"
    mock$Strain=gsub("[[:punct:]].+","",mock$Strain)
    abs=c("IgM","TotalG","IgG1","IgG2a","IgG2b","IgG2c","IgG3","IgA","IgG2ac")

    # dat=dat.15
    # abs=c("IgG1" ,  "IgG2ac" , "IgG2b" , "IgG3" ,  "IgM"  ,  "TotalG")
    
    hc=3
    i=1

    dat=mock
    phenotype=dat[c(1:hc,i+hc)]
    phenotype=phenotype[complete.cases(phenotype),]

    g<-ggplot(phenotype, aes_string(x=abs[i],y=paste0("reorder(Strain,",abs[i],")")))
    g+geom_point(aes(colour=Strain))+theme_classic()+
      labs(title=paste(abs[i],sep=" "),y="Strain",x="µg/µL")

    
### genotype vs. phenotype for mock QTL
    mock$ID=as.factor(mock$ID)
    mock.avg=summaryBy(. ~ Strain, data=mock, FUN=mean, na.rm=T)
    colnames(mock.avg)[2:10] = gsub(".mean","",colnames(mock.avg[2:10]))
    # strain.list=c("CC004","CC030","CC032","CC034","CC062","CC068","CC074","CC073")
    # mock.avg[which(mock.avg$Strain %in% strain.list),c(1,3)]
    
    #averages per RIX
    ms=merge(mock.avg,chr.status,all=T)
    ms=ms[c(1,3,12)]
    g<-ggplot(ms, aes(x=status,y=TotalG))
    g+geom_jitter(aes(colour=status),na.rm=T)+theme_classic()+theme(legend.position='none')
    
    #individual mice
    ms=merge(mock,chr.status,all=T)
    g<-ggplot(ms, aes(x=status,y=TotalG))
    g+geom_jitter(aes(colour=status),na.rm=T)+theme_classic()+theme(legend.position='none')

    
#####compare IAV-induced to baseline allele
    chr.region="chr18_70-82"
    chr.status=read.csv(file=paste("~/Dropbox/Heise/U19-Ab/antibody/data_bin/mock/",chr.region,
                                   "_scores.csv",sep=""))
    colnames(chr.status)[1]="strain"
    
    day=15
    dat=get(paste0("dat.",day))
    
        dat$RIX=as.character(dat$RIX)
        
        RIX_sep <- data.frame(do.call("rbind", strsplit(dat$RIX,"x")))
        colnames(RIX_sep)[1:2]=c("dam","sire")
            dat=cbind(RIX_sep,dat)
        
        qtl=chr.status
        qtl=qtl[c(1,3)]
        qtl[2]=round(qtl[2],digit=2)
      
        colnames(qtl)[1]="dam"
        dat=merge(dat,qtl)
        colnames(dat)[length(dat)]="dam.qtl"
        
        colnames(qtl)[1]="sire"
        dat=merge(dat,qtl)
        colnames(dat)[length(dat)]="sire.qtl"
        
        dat$sum.qtl=rowSums(dat[c("dam.qtl","sire.qtl")])
        dat$sum.qtl=as.factor(dat$sum.qtl)
    
    inheritance.model="sum.qtl"
    response.var="TotalG"
        
    g=ggplot(dat, aes_string(inheritance.model, response.var))
    g+labs(x=paste("Allele Score at Baseline Ab Total IgG QTL under",inheritance.model,"model",sep=" "),y=paste("IAV-Induced",response.var,"on day",day,sep=" ")) +
    geom_jitter(width=0.3,aes(color=RIX),na.rm=T)+
      theme_minimal()+theme(legend.position='none')
    
    # g+stat_boxplot(aes(middle=mean(data), color=sum.qtl)) 
    
    
#Compare based on actual ab concentration, not an allele
day=15
dat=get(paste0("dat.",day))

### code with  baseline ab (e.g. for graphing)
dat$RIX=as.character(dat$RIX)
    
RIX_sep <- data.frame(do.call("rbind", strsplit(dat$RIX,"x")))
colnames(RIX_sep)[1:2]=c("dam","sire")
    
dat=cbind(RIX_sep,dat)

ug=mock.avg[c(1,3)]
colnames(ug)[2]=paste(colnames(ug)[2],".RI",sep="")

colnames(ug)[1]="sire"
dat=merge(dat,ug,by="sire")
colnames(dat)[length(dat)]="sire.ug"

colnames(ug)[1]="dam"
dat=merge(dat,ug,by="dam")
colnames(dat)[length(dat)]="dam.ug"

dat$mean.ab=rowMeans(dat[c("dam.ug","sire.ug")])

inheritance.model="mean.ab"
response.var="TotalG"
    
g=ggplot(dat, aes_string(inheritance.model, response.var))
g+labs(x="Sum of Parent RI Total IgG at Baseline (µg/mL)",y=paste("IAV-Induced",response.var,"on day",day,sep=" ")) +
  geom_jitter(width=0.1,aes(colour = RIX),na.rm=T)+theme_minimal()+theme(legend.position='none')

lm1=lm(TotalG~mean.ab,data=dat)
summary(lm1)
    