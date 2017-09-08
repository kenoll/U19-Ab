# install.packages("data.table")
library("data.table")

setwd("~/Dropbox/Heise")

select_gene_ids <- c("ENSMUSG00000070390", "ENSMUSG00000069830", "ENSMUSG00000000386",
                     "ENSMUSG00000032691","ENSMUSG00000040296","ENSMUSG00000027514")
select_gene_labels = c("Nlrp1b","Nlrp1a","Mx1","Nlrp3","RIG-I","Zbp1")

scaled <- as.data.frame(fread("GSE52405_scaled_counts_plm.csv", header=TRUE))
rownames(scaled) <- as.character(scaled$V1)
scaled <- scaled[,-1]

founder.colors <- c("#F0F000", "#808080", "#F08080", "#1010F0", "#00A0F0", "#00A000", "#F00000", "#9000E0")
founder.names <- c("AJ", "C57BL6J", "129S1", "NOD", "NZO", "CAST", "PWK", "WSB")

select_scaled <- scaled[,as.character(colnames(scaled)) %in% c(as.character(select_gene_ids), "Strain", "Trt", "Day", "SampleID")]
select_scaled <- subset(select_scaled, !(Trt=="MA15"))
select_scaled$Strain <- factor(select_scaled$Strain, levels=founder.names)
select_scaled$color <- founder.colors[match(select_scaled$Strain, founder.names)]
select_scaled$DayTrtStrain <- factor(paste(select_scaled$Day, select_scaled$Trt, select_scaled$Strain,sep="."))#, levels=reorder)
ids <- select_gene_ids[select_gene_ids %in% colnames(scaled)]

pdf("~/Desktop/cc-founder-expression-contrast.pdf", width=11, height=8.5)
par(mfrow=c(2,3))
par(mar=c(10,3,3,1))
# select_scaled_contrast <- droplevels(subset(select_scaled, Strain %in% c( "AJ", "C57BL6J", "129S1", "PWK")))
for(i in 1:length(ids)){
# try(main <- gene_list_Chr2_all[gene_list_Chr2_all$ensembl_gene_id==ids[i],]$mgi_symbol)
  boxplot(select_scaled[,i]~DayTrtStrain, data=select_scaled, las=2, main=select_gene_labels[i],
          col=founder.colors[c(3,1,2,6,4,5,7,8)])
  abline(v=0.5+8*c(1:3))
}
dev.off()

#another example
# 
# par(mfrow=c(1,1))
# scaled$DayTrtStrain <- factor(paste(scaled$Day, scaled$Trt, scaled$Strain,sep="."))
# scaled <- droplevels(subset(scaled, !(Trt=="MA15")))
# main <- "Eif3m"
# boxplot(ENSMUSG00000070390~DayTrtStrain, data=scaled, las=2, main=main,
#         col=founder.colors[c(3,1,2,6,4,5,7,8)])
# abline(v=0.5+8*c(1:3))
