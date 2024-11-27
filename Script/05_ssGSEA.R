rm(list = ls())
#引用包
library(reshape2)
library(ggpubr)
library(limma)
library(GSEABase)
library(GSVA)

gmtFile="immune.gmt"          #免疫数据集文件
clusterFile="Cluster.txt"       #分型输入文件

#读取表达输入文件,并对输入文件整理
load("cleandata.RDATA")

exp=brca.tcga$mRNA.expr
data=exp

#读取基因集文件
geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())

#ssGSEA分析
data=as.matrix(data)
ssgseaScore=gsva(data, geneSets, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
#对ssGSEA打分进行矫正
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaScore=normalize(ssgseaScore)
#输出ssGSEA打分结果
ssgseaOut=rbind(id=colnames(ssgseaScore), ssgseaScore)
write.table(ssgseaOut,file="ssGSEA.result.txt",sep="\t",quote=F,col.names=F)

#读取分型文件
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

#数据合并
ssgseaScore=t(ssgseaScore)
library(tidyverse)


sameSample=intersect(row.names(ssgseaScore), row.names(cluster))



ssgseaScore=ssgseaScore[sameSample,,drop=F]
cluster=cluster[sameSample,,drop=F]
scoreCluster=cbind(ssgseaScore, cluster)

#把数据转换成ggplot2输入文件
data=melt(scoreCluster, id.vars=c("cluster"))
colnames(data)=c("cluster", "Immune", "Fraction")

#绘制箱线图
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data[,"cluster"])))]
library(tidyverse)
data$Immune=str_replace_all(data$Immune,"na","")
p=ggboxplot(data, x="Immune", y="Fraction", color="cluster", 
            ylab="Immune infiltration",
            xlab="",
            legend.title="cluster",
            palette=bioCol)
p=p+rotate_x_text(50)
pdf(file="boxplot.pdf", width=10, height=6.5)                          #输出图片文件
p+stat_compare_means(aes(group=cluster),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
dev.off()



