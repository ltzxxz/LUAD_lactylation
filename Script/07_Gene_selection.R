#加载包
rm(list = ls())
library(tidyverse)
library(ggpubr)
library(survminer)
library(survival)
library(survivalROC)
library(pheatmap)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
library(ReactomePA)
library(igraph)
library(ggraph)
library(ggradar)
library(tuneR)
library(limma)
library(stringr)
library(RColorBrewer)
library(forestplot)
library(fmsb)
library(circlize)
library(ggsci)
library(parallel)
library(maftools)
library(circlize)
library(ggsci)
library(parallel)
library(patchwork)

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
rm(list = ls())
load("cleandata.RDATA")
load("GSE13213.RDATA")
load("GSE31210.RDATA")


Index<-c("TCGA","GSE13213","GSE31210")

genename=read.table("markergene.txt",header = T,sep = "\t")
genename=genename$probe
outTab=brca.tcga$mRNA.expr
outTab=outTab[genename,]

aaa=as.data.frame(t(outTab))
aaa$ID=rownames(aaa)
aaa=dplyr::select(aaa,ID,everything())

clin=brca.tcga$clin.info
clin$ID=rownames(clin)
clin=dplyr::select(clin,ID,everything())
# aaa$ID=str_replace_all(aaa$ID,"TCGA_","")
# aaa$ID=str_replace_all(aaa$ID,"GSE13213_","")
OSdata=inner_join(clin,aaa,by="ID")
colnames(OSdata)[2:3]=c("OS.time","OS")
#OS


library(survival)   
library(tidyverse)
OSdata$OS.time=OSdata$OS.time/365  
OSdata1=OSdata

rt<-dplyr::select(OSdata,ID,OS,OS.time,genename)
rt=dplyr::filter(rt,OS.time > 0)
 

outTab=data.frame()
for(gene in    colnames(rt)[4:ncol(rt)]       ){
  print(which(gene== colnames(rt)[4:ncol(rt)]  ))
  cox=coxph(Surv(OS.time, OS) ~ rt[,gene], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(gene=gene,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxP) )
}
outTab$HR=as.numeric(outTab$HR)
outTab$HR.95L=as.numeric(outTab$HR.95L)
outTab$HR.95H=as.numeric(outTab$HR.95H)
outTab$pvalue=as.numeric(outTab$pvalue)
write.table(outTab,file=paste0("TCGA","单因素回归分析结果.txt"),sep="\t",row.names=F,quote=F)
badgene=dplyr::filter(outTab,pvalue < 0.05,HR>1)
assign(paste0("TCGA",".badgene"),badgene$gene)
goodgene=dplyr::filter(outTab,pvalue < 0.05,HR<1)
assign(paste0("TCGA",".goodgene"),goodgene$gene)

#######################################################################
#######################################################################
#######################################################################

genename=read.table("markergene.txt",header = T,sep = "\t")
genename=genename$probe
outTab=GSE13213$mRNA.expr

genename=intersect(genename,rownames(outTab))
outTab=outTab[genename,]

aaa=as.data.frame(t(outTab))
aaa$ID=rownames(aaa)
aaa=dplyr::select(aaa,ID,everything())

clin=GSE13213$clin.info
clin$ID=rownames(clin)
clin=dplyr::select(clin,ID,everything())
# aaa$ID=str_replace_all(aaa$ID,"TCGA_","")
# aaa$ID=str_replace_all(aaa$ID,"GSE13213_","")
OSdata=inner_join(clin,aaa,by="ID")
colnames(OSdata)[2:3]=c("OS.time","OS")
#OS


library(survival)   
library(tidyverse)
OSdata$OS.time=OSdata$OS.time/365  
OSdata2=OSdata

rt<-dplyr::select(OSdata,ID,OS,OS.time,genename)
rt=dplyr::filter(rt,OS.time > 0)


outTab=data.frame()
for(gene in    colnames(rt)[4:ncol(rt)]       ){
  print(which(gene== colnames(rt)[4:ncol(rt)]  ))
  cox=coxph(Surv(OS.time, OS) ~ rt[,gene], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(gene=gene,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxP) )
}
outTab$HR=as.numeric(outTab$HR)
outTab$HR.95L=as.numeric(outTab$HR.95L)
outTab$HR.95H=as.numeric(outTab$HR.95H)
outTab$pvalue=as.numeric(outTab$pvalue)
write.table(outTab,file=paste0("GSE13213","单因素回归分析结果.txt"),sep="\t",row.names=F,quote=F)
badgene=dplyr::filter(outTab,pvalue < 0.05,HR>1)
assign(paste0("GSE13213",".badgene"),badgene$gene)
goodgene=dplyr::filter(outTab,pvalue < 0.05,HR<1)
assign(paste0("GSE13213",".goodgene"),goodgene$gene)


#######################################################################
#######################################################################
#######################################################################
genename=read.table("markergene.txt",header = T,sep = "\t")
genename=genename$probe
outTab=GSE31210$mRNA.expr

genename=intersect(genename,rownames(outTab))
outTab=outTab[genename,]

aaa=as.data.frame(t(outTab))
aaa$ID=rownames(aaa)
aaa=dplyr::select(aaa,ID,everything())

clin=GSE31210$clin.info
clin$ID=rownames(clin)
clin=dplyr::select(clin,ID,everything())
# aaa$ID=str_replace_all(aaa$ID,"TCGA_","")
# aaa$ID=str_replace_all(aaa$ID,"GSE13213_","")
OSdata=inner_join(clin,aaa,by="ID")
colnames(OSdata)[2:3]=c("OS.time","OS")
#OS


library(survival)   
library(tidyverse)
OSdata$OS.time=OSdata$OS.time/365  
OSdata3=OSdata

rt<-dplyr::select(OSdata,ID,OS,OS.time,genename)
rt=dplyr::filter(rt,OS.time > 0)


outTab=data.frame()
for(gene in    colnames(rt)[4:ncol(rt)]       ){
  print(which(gene== colnames(rt)[4:ncol(rt)]  ))
  cox=coxph(Surv(OS.time, OS) ~ rt[,gene], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(gene=gene,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxP) )
}
outTab$HR=as.numeric(outTab$HR)
outTab$HR.95L=as.numeric(outTab$HR.95L)
outTab$HR.95H=as.numeric(outTab$HR.95H)
outTab$pvalue=as.numeric(outTab$pvalue)
write.table(outTab,file=paste0("GSE31210","单因素回归分析结果.txt"),sep="\t",row.names=F,quote=F)
badgene=dplyr::filter(outTab,pvalue < 0.05,HR>1)
assign(paste0("GSE31210",".badgene"),badgene$gene)
goodgene=dplyr::filter(outTab,pvalue < 0.05,HR<1)
assign(paste0("GSE31210",".goodgene"),goodgene$gene)



# 创建向量列表
#vecs <- list(CGGA_301.badgene, CGGA_325.badgene, CGGA_693.badgene,GSE13041.badgene,Rembrandt__475.badgene,TCGA_702.badgene)
vecs <- list(TCGA.badgene, GSE13213.badgene,GSE31210.badgene)
#选择至少在4个cohort中出现的交集基因
# 将所有向量放入列表

# 初始化结果为空向量
result <- c()
aaa=c(TCGA.badgene, GSE13213.badgene,GSE31210.badgene)
aaa=unique(aaa)
# 循环遍历向量的每个元素
for (elem in aaa) {
  # 统计该元素在多少个向量中出现
  count <- sum(sapply(vecs, function(x) elem %in% x))
  
  # 如果在至少4个向量中出现，则加入结果向量
  if (count >= 3) {
    result <- c(result, elem)
  }
}

# 打印结果
print(result)
badgene=result


# 创建向量列表
vecs <- list(TCGA.goodgene, GSE13213.goodgene,GSE31210.goodgene)

#选择至少在4个cohort中出现的交集基因
# 将所有向量放入列表

# 初始化结果为空向量
result <- c()
aaa=c(TCGA.goodgene, GSE13213.goodgene,GSE31210.goodgene)
aaa=unique(aaa)
# 循环遍历向量的每个元素
for (elem in aaa) {
  # 统计该元素在多少个向量中出现
  count <- sum(sapply(vecs, function(x) elem %in% x))
  
  # 如果在至少4个向量中出现，则加入结果向量
  if (count >= 3) {
    result <- c(result, elem)
  }
}

# 打印结果
print(result)
goodgene=result


sss=Index[1]
rt=read.table(paste0(sss,"单因素回归分析结果.txt"),header=T,sep="\t",row.names=1,check.names=F)
rt=rt[c(badgene,goodgene),]
rt=arrange(rt,desc(HR))
#把含有0或者NA的行全都去掉
rt[rt == 0] <- NA
rt<-na.omit(rt)
options(forestplot_new_page = FALSE)
clrs <- fpColors(box="red",line="darkblue", summary="royalblue") 
data=as.matrix(rt)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )   
pdf(file=paste0(sss,"_OS_forest.pdf"),
    width = 8,            
    height = 12,           
)

forestplot(tabletext, 
           zero = 1,
           lwd.zero = 2,
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(50, "mm"),
           xlog=T,
           lwd.ci=2,
           boxsize=0.4,
           xlab="Hazard ratio"
)
dev.off()

OSdata1=dplyr::filter(OSdata1,OS.time > 0)

OSdata1=OSdata1[,c(colnames(OSdata)[1:10],badgene,goodgene)]
OSdata1$cohort="TCGA"
OSdata2=OSdata2[,c(colnames(OSdata)[1:10],badgene,goodgene)]
OSdata2$cohort="GSE13213"
OSdata3=OSdata3[,c(colnames(OSdata)[1:10],badgene,goodgene)]
OSdata3$cohort="GSE31210"
OSdata=rbind(OSdata1,OSdata2,OSdata3)
OSdata=dplyr::select(OSdata,ID,OS.time,OS,cohort,everything())

save(goodgene,badgene,OSdata,file =  "筛选基因表达临床.RDATA")
