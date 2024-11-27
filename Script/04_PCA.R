rm(list = ls())
#引用包
library(limma)
library(ggplot2)

clusterFile="Cluster.txt"     #分型文件

#读取输入文件,并对输入文件进行整理
load("cleandata.RDATA")
rt=brca.tcga$mRNA.expr
rt=as.matrix(rt)
exp=rt
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
# data=data[rowMeans(data)>0,]
data=t(data)

#PCA分析
data.pca=prcomp(data, scale. = F)
pcaPredict=predict(data.pca)
write.table(pcaPredict, file="newTab.xls", quote=F, sep="\t")

#读取分型文件
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
cluster=as.vector(cluster[,1])

#设置颜色
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
CluCol=bioCol[1:length(levels(factor(cluster)))]

#可视化
PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], cluster=cluster)
PCA.mean=aggregate(PCA[,1:2], list(cluster=PCA$cluster), mean)
pdf(file="PCA.pdf", height=5, width=6.5)
ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = cluster)) +
  scale_colour_manual(name="cluster", values =CluCol)+
  theme_bw()+
  theme(plot.margin=unit(rep(1.5,4),'lines'))+
  annotate("text",x=PCA.mean$PC1, y=PCA.mean$PC2, label=PCA.mean$cluster, cex=7)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()


