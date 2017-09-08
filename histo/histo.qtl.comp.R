histo=read.csv("histo/flu-histo-map.csv")
pheno=read.csv("weight/peak_weight_loss.csv")

rixs=read.csv("qtl_info/chr1_38-41_U19RIXhaplos.csv")
score=read.csv("qtl_info/chr1_38-41_scores.csv")

histo=histo[c(3,6:20)]

colnames(histo)[1]<-"RIX"
histo$RIX=as.character(histo$RIX)
histo=alias.to.line(histo)
histo=add.stat.gene(histo,score,"chr1")
histo=add.stat.mx1(histo)

histo.mx1=subset(histo,histo$add_mx1==0)
table(histo.mx1$add_chr1)

histo.chr1=subset(histo.mx1,histo.mx1$add_chr1==3 | histo.mx1$add_chr1==1)
histo.chr1=histo.chr1[c(1:18,21,25)]

abs=colnames(histo.chr1[7:18])

library(lme4)

for (i in 1:length(abs)) {
  
test=abs[i]

model_base<-lmer(get(test)~(1|RIX),data=histo.chr1)
model_full<-lmer(get(test)~(1|RIX)+add_chr1, data=histo.chr1)
anv=anova(model_base, model_full)
print(anv)
}
