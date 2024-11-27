rm(list = ls())
#加载包
library(tidyverse)
library(maftools)
load("RS_Risk_clin.RDATA")
drawdata<-OSdata

laml.maf=data.table::fread(paste0("TCGA.","LUAD",".varscan.gz"),data.table = F)
laml.maf$Tumor_Sample_Barcode<-substr(laml.maf$Tumor_Sample_Barcode,1,12)
laml.maf$Tumor_Sample_Barcode<-str_replace_all(laml.maf$Tumor_Sample_Barcode,"-",".")

drawdata=dplyr::filter(drawdata,cohort == "TCGA")

drawdata$ID=substr(drawdata$ID,1,12)
drawdata$ID=str_replace_all(drawdata$ID,"_",".")



newdata=drawdata

value<-newdata$ID[which(newdata$ID %in% laml.maf$Tumor_Sample_Barcode)]
newdata<-newdata%>%
  dplyr::filter(ID %in%value )


Highexpr<-newdata$ID[which(newdata[,"Risk_group"]=="High_risk")]
Lowexpr<-newdata$ID[which(newdata[,"Risk_group"]=="Low_risk")]

High.laml.maf<-laml.maf%>%
  dplyr::filter(Tumor_Sample_Barcode %in% Highexpr )
Low.laml.maf<-laml.maf%>%
  dplyr::filter(Tumor_Sample_Barcode %in% Lowexpr )

High.laml = read.maf(maf = High.laml.maf)
Low.laml = read.maf(maf = Low.laml.maf)

vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)

pdf(file=paste0("高score组突变瀑布图.pdf"),width=8,height=6)
oncoplot(maf = High.laml, colors = vc_cols,top = 10)
dev.off()
pdf(file=paste0("低score组突变瀑布图.pdf"),width=8,height=6)
oncoplot(maf = Low.laml, colors = vc_cols, top = 10)
dev.off()
pt.vs.rt <- mafCompare(m1 = High.laml, m2 = Low.laml, m1Name = paste0("High_risk"), m2Name =  paste0("Low_risk"), minMut = 0)
write.csv(pt.vs.rt$results,file =paste0("高低score组差异突变基因对比.csv"))
#选取前12个基因展示
pt.vs.rt$results<-pt.vs.rt$results[1:12,]
pdf(file=paste0("高低score组差异突变基因对比.pdf"),width=8,height=10)
maftools::forestPlot(mafCompareRes = pt.vs.rt, pVal = 1, geneFontSize = 0.8)
dev.off()

  

