rm(list = ls())
library(RColorBrewer)
library(circlize)
library(gplots)
library(viridis)
library(oompaBase)
library(tidyverse)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
  outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}
load("cleandata.RDATA")
load("RS_Risk_clin.RDATA")
dataa=as.data.frame(brca.tcga$mRNA.expr)

dataa=dataa[,OSdata$ID]
aaa=read.csv("Cytokine.csv",header = T,sep = ",")
genename1=intersect(aaa$ID,rownames(dataa))
dataa=dataa[genename1,]
aaa=dplyr::filter(aaa,ID %in% genename1)


risk=OSdata
rownames(risk)=risk$ID
sameid=intersect(colnames(dataa),rownames(risk))
dataa=dataa[,sameid]
risk=risk[sameid,]
hmdat=as.data.frame(t(dataa))
#分类信息
type=data.frame(Type =aaa$Type)
row.names(type)=aaa$ID

Typeid <- type$Type
colnames(risk)
# 创建注释
# 列注释，位于热图顶端
annCol <- data.frame(score = scale(risk$RS),
                     RiskType = risk$Risk_group,
                     # 以上是risk score和risk type两种注释，可以按照这样的格式继续添加更多种类的注释信息，记得在下面的annColors里设置颜色
                     row.names = rownames(risk),
                     stringsAsFactors = F)


#基因差异分析
bbb=annCol[rownames(hmdat),]
sampleType=bbb$RiskType
sigVec=c()
for(i in colnames(hmdat)){
  test=wilcox.test(hmdat[,i] ~ sampleType)
  pvalue=test$p.value
  Sig=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
  sigVec=c(sigVec, paste0(i, Sig))
}

colnames(hmdat)=sigVec
# 数据标准化

# 用pheatmap画图
library(pheatmap)

# 定义颜色
Type.col <- brewer.pal(n = length(unique(Typeid)),name = "Paired")


# 行注释，位于热图左侧
annRow <- data.frame(Type = factor(Typeid,levels = unique(Typeid)),
                     row.names = colnames(hmdat),
                     stringsAsFactors = F)


unique(annRow$Type)
# 为各注释信息设置颜色
annColors <- list(Type = c("Chemokines and receptors" = Type.col[1], #行注释的颜色
                           "Interieukines and receptors" = Type.col[2],
                           "Interferons and receptors" = Type.col[3],
                           "Other cytokines" = Type.col[4]),
                  # 下面是列注释的颜色，可依此设置更多注释的颜色
                  "RS" = greenred(64), 
                  "RiskType" = c("High_risk" = "red","Low_risk" = "blue"))



# 样本按risk score排序
samorder <- rownames(risk[order(risk$RS),])



indata <- t(hmdat)
indata <- indata[,colSums(indata) > 0] # 确保没有富集全为0的细胞
plotdata <- standarize.fun(indata,halfwidth = 2)

# pheatmap绘图
pheatmap::pheatmap(mat = as.matrix(plotdata[,samorder]), # 标准化后的数值矩阵
                   border_color = NA, # 无边框色
                   color = bluered(64), # 热图颜色为红蓝
                   cluster_rows = F, # 行不聚类
                   cluster_cols = F, # 列不聚类
                   show_rownames = T, # 显示行名
                   show_colnames = F, # 不显示列名
                   annotation_col = annCol[samorder,,drop = F], # 列注释
                   annotation_row = annRow, # 行注释
                   annotation_colors = annColors, # 注释颜色
                   gaps_col = table(annCol$RiskType)[2], # 列分割
                   gaps_row = cumsum(table(annRow$Type)), # 行分割
                   cellwidth = 1, # 元素宽度
                   cellheight = 10, # 元素高度
                   filename = "immune heatmap by pheatmap_score.pdf")




hmdat$ID=rownames(hmdat)
annCol$ID=rownames(annCol)
annCol=dplyr::select(annCol,ID,everything())
newdata=inner_join(annCol,hmdat,by="ID")


inputtemp<-newdata

colnames(inputtemp) <- gsub("\\**", "", colnames(inputtemp))
colnames(inputtemp) <- gsub("\\*", "", colnames(inputtemp))

gene="score"
value=colnames(inputtemp)[4:106]
length(value)
for (b in 1:length(value)) {
  test <- cor.test(inputtemp[,gene],inputtemp[,value[b]],exact=FALSE)
  paste(paste0("n = ",length(inputtemp[,gene])),
        paste0("r = ",round(test$estimate,2),"(",'pearson',")"),
        paste0("p.value= ",round(test$p.value,4)),
        sep = ", ")
  p<-ggplot(inputtemp,aes(get(gene),get(value[b])))+
    geom_point(col="#984ea3")+
    geom_smooth(method=lm, se=T,na.rm=T, fullrange=T,size=2,col="#fdc086")+
    geom_rug(col="#7fc97f")+
    theme_minimal()+
    xlab(paste("Riskscore",sep = ' '    ))+
    ylab(paste(value[b],sep = ' '    ))+
    ## 依靠函数来生成title
    labs(title = paste(
      paste0("n = ",length(inputtemp[,gene])),
      paste0("r = ",round(test$estimate,2),"(",'pearson',")"),
      paste0("p.value= ",round(test$p.value,4)),
      sep = ", "))+
    theme(plot.title = element_text(hjust = 0.5),
          plot.margin = margin(1, 1, 1, 1, "cm"))
  ggsave(p,filename = paste(gene,value[b],'相关性点图.pdf',sep=' '),width = 4.94,height = 4.72)
}



















