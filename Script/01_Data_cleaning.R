rm(list = ls())
library(tidyverse)
library(data.table)
library(maftools)
library("MOVICS")

#mRNA表达谱
mRNA.expr=read.table("TCGA.txt",header = T,sep = "\t",row.names = 1)
str(mRNA.expr)
#lncRNA表达谱
lcRNA.expr=read.table("TCGA-LUAD.LNCRNA.TXT",header = T,sep = "\t",row.names = 1)
str(lcRNA.expr)

#突变数据
laml.maf <- data.table::fread("TCGA.LUAD.varscan.gz",data.table = F)
laml.maf$ID=substr(laml.maf$Tumor_Sample_Barcode,1,16)
allgene=unique(laml.maf$Hugo_Symbol)

# 指定数据框的行数和列数
rows <- length(allgene)
cols <- length(unique(laml.maf$ID))

# 创建一个全为NA的矩阵
mat <- matrix(0, nrow = rows, ncol = cols)

# 将矩阵转换为数据框
mut.status <- as.data.frame(mat)
colnames(mut.status)=unique(laml.maf$ID)
rownames(mut.status)=allgene
dim(mut.status)

for (i in 1:length(unique(laml.maf$ID))) {
  print(i)
  sampleid=unique(laml.maf$ID)[i]
  adata=dplyr::filter(laml.maf,ID  %in% sampleid)
  mut.status[adata$Hugo_Symbol,sampleid]=1
}
str(mut.status)
#甲基化数据
meth.beta=read.table("data_methylation_hm27_hm450_merged.txt",header = T,sep = "\t",row.names = 1)
meth.beta=na.omit(meth.beta)
meth.beta=as.matrix(meth.beta)
rownames(meth.beta)=meth.beta[,"NAME"]
meth.beta=meth.beta[,-c(1,2,3)]

library(limma)
meth.beta=avereps(meth.beta)
meth.beta=as.data.frame(meth.beta)
meth.beta1=meth.beta
meth.beta <- as.data.frame(lapply(meth.beta, as.numeric))
rownames(meth.beta)=rownames(meth.beta1)
rownames(meth.beta)=str_replace_all(rownames(meth.beta),"-",".")
rownames(meth.beta)=str_replace_all(rownames(meth.beta),";",".")




#临床资料
clin.info=read.table("clinicaldata.txt",header = T,sep = "\t")
clinid=intersect(colnames(mRNA.expr),clin.info$ID  )
clin.info=dplyr::filter(clin.info, ID %in% clinid)

#处理ID
colnames(mRNA.expr)=substr(colnames(mRNA.expr),1,15)
colnames(lcRNA.expr)=substr(colnames(lcRNA.expr),1,15)
colnames(mut.status)=substr(colnames(mut.status),1,15)
clin.info$ID=substr(clin.info$ID,1,15)


colnames(mRNA.expr)=str_replace_all(colnames(mRNA.expr),"-",".")
colnames(lcRNA.expr)=str_replace_all(colnames(lcRNA.expr),"-",".")
colnames(mut.status)=str_replace_all(colnames(mut.status),"-",".")
clin.info$ID=str_replace_all(clin.info$ID,"-",".")



aaa=intersect(colnames(mRNA.expr),colnames(lcRNA.expr))
bbb=intersect(aaa,colnames(mut.status))
ccc=intersect(bbb,colnames(meth.beta))
ddd=intersect(ccc,clin.info$ID)


sameid=ddd
mRNA.expr=mRNA.expr[,sameid]
lcRNA.expr=lcRNA.expr[,sameid]
mut.status=mut.status[,sameid]
meth.beta=meth.beta[,sameid]
clin.info=dplyr::filter(clin.info,ID %in% sameid)
clin.info=distinct(clin.info,ID,.keep_all = T)
rownames(clin.info)=clin.info$ID
clin.info=clin.info[,-1]
clin.info=clin.info[sameid,]
brca.tcga=list(mRNA.expr=mRNA.expr,lcRNA.expr=lcRNA.expr,meth.beta=meth.beta,mut.status=mut.status,clin.info=clin.info)

save(brca.tcga,file =  "cleandata.RDATA")


rm(list = ls())
library(tidyverse)
library(data.table)
library(maftools)
library("MOVICS")

#mRNA表达谱
mRNA.expr=read.table("GSE31210.txt",header = T,sep = "\t",row.names = 1)
str(mRNA.expr)


#临床资料
clin.info=read.table("GSE31210.clin.txt",header = T,sep = "\t")
clinid=intersect(colnames(mRNA.expr),clin.info$ID  )
clin.info=dplyr::filter(clin.info, ID %in% clinid)


ddd=intersect(colnames(mRNA.expr),clin.info$ID)



sameid=ddd
mRNA.expr=mRNA.expr[,sameid]
clin.info=dplyr::filter(clin.info,ID %in% sameid)
clin.info=distinct(clin.info,ID,.keep_all = T)
rownames(clin.info)=clin.info$ID
clin.info=clin.info[,-1]
clin.info=clin.info[sameid,]
GSE31210=list(mRNA.expr=mRNA.expr,clin.info=clin.info)
save(GSE31210,file =  "GSE31210.RDATA")





rm(list = ls())
library(tidyverse)
library(data.table)
library(maftools)
library("MOVICS")

#mRNA表达谱
mRNA.expr=read.table("GSE13213.txt",header = T,sep = "\t",row.names = 1)
str(mRNA.expr)


#临床资料
clin.info=read.table("GSE13213.clin.txt",header = T,sep = "\t")
clinid=intersect(colnames(mRNA.expr),clin.info$ID  )
clin.info=dplyr::filter(clin.info, ID %in% clinid)


ddd=intersect(colnames(mRNA.expr),clin.info$ID)



sameid=ddd
mRNA.expr=mRNA.expr[,sameid]
clin.info=dplyr::filter(clin.info,ID %in% sameid)
clin.info=distinct(clin.info,ID,.keep_all = T)
rownames(clin.info)=clin.info$ID
clin.info=clin.info[,-1]
clin.info=clin.info[sameid,]
GSE13213=list(mRNA.expr=mRNA.expr,clin.info=clin.info)
save(GSE13213,file =  "GSE13213.RDATA")
