setwd("~/Dropbox/Heise/U19-Ab")
library(reshape2)

sars.ab=read.csv("sars_ab/data_bin/SARS_halfmax_lastpos.csv")
flu.ab=read.csv("antibody/data_bin/FLU_halfmax_lastpos.csv")
abs=c("IgG1"  , "IgG2ac" , "IgG2b" ,  "IgG3"  ,  "IgM"  ,  "TotalG")


flu.ab$day=paste0("D",flu.ab$day)

sars.ab=sars.ab[c(1,3:9)]
flu.ab=flu.ab[c(1:2,4:5,7,9,11:12)]

sars.half=dcast(sars.ab, RIX + ID + day + virus + antigen ~ isotype,mean,value.var='halfmax')
flu.half=dcast(flu.ab, RIX + ID + day + virus + antigen ~ isotype,mean,value.var='halfmax')

half.dat=rbind(sars.half,flu.half)

sars.last=dcast(sars.ab, RIX + ID + day + virus + antigen ~ isotype,mean,value.var='last_positive')
flu.last=dcast(flu.ab, RIX + ID + day + virus + antigen ~ isotype,mean,value.var='last_positive')

last.dat=rbind(sars.last,flu.last)

last.avg=summaryBy(IgG1+IgG2ac+IgG2b+IgG3+IgM+TotalG~RIX+day+virus+antigen,data=last.dat,FUN=mean, na.rm=T)
colnames(last.avg)[5:10] = gsub(".mean","",colnames(last.avg[5:10]))


####

ab.dat=rbind(sars.ab,flu.ab)

plot.data=subset(ab.dat,ab.dat$day=="D7")

strain.abs=summaryBy(halfmax+last_positive~RIX+day+isotype+virus+antigen,data=ab.dat,FUN=mean, na.rm=T)
colnames(strain.abs)[(length(strain.abs)-1):length(strain.abs)]=c("halfmax","last_positive")


###
colnames(sars.ab)[7:8]=c("sars.half","sars.last")
colnames(flu.ab)[7:8]=c("flu.half","flu.last")

###

model_base<-lmer(get(test)~(1|RIX),data=histo.chr1)
model_full<-lmer(get(test)~(1|RIX)+add_chr1, data=histo.chr1)
anv=anova(model_base, model_full)
print(anv)



head(sars.10)
chr11=read.csv("antibody/phenotypes/d10/chr11_70.5-72_scores.csv")

sars.10 = sars.10 %>% add.stat.gene(chr11,"chr11")

# 
# ggplot(plot.data, aes(x=reorder(RIX,value),y=value,group=virus,color=antigen))+
#   labs(x="RIX",y=paste0(isotype))+
#   theme_minimal()+
#   geom_point(na.rm=T)+
#   theme(axis.text.x=element_text(angle=60,vjust=1,hjust=1))+scale_x_discrete(drop=FALSE)
# 







# ##############
# 
# sars.7=read.csv("sars_ab/data_bin/SARS_D7_lastpos.csv")
# sars.10=read.csv("sars_ab/data_bin/SARS_D10_lastpos.csv")
# sars.15=read.csv("sars_ab/data_bin/SARS_D15_lastpos.csv")
# 
# flu.7=read.csv("antibody/data_bin/d7/dat.7.lastpos.csv")
# flu.10=read.csv("antibody/data_bin/d10/dat.10.lastpos.csv")
# flu.15=read.csv("antibody/data_bin/d15/dat.15.lastpos.csv")
# 
# day=15
# 
# df.sars=get(paste0("sars.",day))
# df.flu=get(paste0("flu.",day))
# 
# df.sars=melt(df.sars,id.vars=colnames(df.sars)[1:6])
# df.flu=melt(df.flu,id.vars=colnames(df.flu)[1:7])
# 
# df.full=rbind(df.sars[c(1,3:8)],df.flu[c(1:5,8:9)])
# df.full$day=day
# 
# isotype="TotalG"
# 
# plot.data=subset(df.full,df.full$variable==isotype)
# 
# ggplot(plot.data, aes(x=reorder(RIX,value),y=value,group=virus,color=antigen))+
#   labs(x="RIX",y=paste0(isotype))+
#   theme_minimal()+
#   geom_point(na.rm=T)+
#   theme(axis.text.x=element_text(angle=60,vjust=1,hjust=1))+scale_x_discrete(drop=FALSE)

