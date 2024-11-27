rm(list = ls())


library(RColorBrewer)
library(circlize)
library(gplots)
library(viridis)
library(oompaBase)
library(tidyverse)
library(pheatmap)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ReactomePA)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor


load("cleandata.RDATA")
outTab=brca.tcga$mRNA.expr
load("RS_Risk_clin.RDATA")
clin=OSdata

outTab=as.data.frame(outTab)

outTab=outTab[,clin$ID]


outTab=as.data.frame(t(outTab))

outTab$ID=rownames(outTab)
clin=clin[,c("ID","RS")]
colnames(clin)[2]="Riskscore"

FF=inner_join(clin,outTab,by="ID")

gene="Riskscore"

y <- as.numeric(FF[,gene])#开始相关性分析

colnames <- colnames(FF)[3:length(FF)]
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(FF[,i+2]),y,type="pearson")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")
write.csv (cor_data_df, file =paste( gene,'批量相关性分析.csv',sep=' '), row.names =FALSE)#将文件导出
#相关性分析结果排序选取P<0.05,top50的做相关性热图，R>0.3的做相关性图
cor_data_df<-na.omit(cor_data_df)
pos<-cor_data_df %>%
  filter(pvalue < 0.05 ) %>%
  arrange(desc(correlation)) %>%
  top_n(50,correlation)
neg<-cor_data_df %>%
  filter(pvalue < 0.05 ) %>%
  arrange(correlation) %>%
  top_n(-49,correlation)
#正相关热图
phemapdata<-FF %>%
  dplyr::select(ID,gene,pos$symbol)
phemapdata<-arrange(phemapdata,desc(phemapdata[,gene]))
row.names(phemapdata)<-phemapdata$ID
phemapdata<-phemapdata[,-1]
phemapdata<-as.data.frame(t(phemapdata))
p<-pheatmap(phemapdata,scale="row",cluster_row=F,cluster_col=F,legend= T,show_colnames = F,color = colorRampPalette(c("blue", "white", "red"))(50))
pdf(file = paste( gene,'正相关性热图.pdf',sep=' '),width = 8,height = 6)
print(p)
dev.off()
#负相关热图
phemapdata<-FF %>%
  dplyr::select(ID,gene,neg$symbol)
phemapdata<-arrange(phemapdata,desc(phemapdata[,gene]))
row.names(phemapdata)<-phemapdata$ID
phemapdata<-phemapdata[,-1]
phemapdata<-as.data.frame(t(phemapdata))
p<-pheatmap(phemapdata,scale="row",cluster_row=F,cluster_col=F,legend= T,show_colnames = F,color = colorRampPalette(c("blue", "white", "red"))(50))
pdf(file = paste( gene,'负相关性热图.pdf',sep=' '),width = 8,height = 6)
print(p)
dev.off()
#基因相关性作图（选取正负相关最显著的各50个基因作图）
#setwd("Y:/lcd代码/TCGA随意做/TCGA任意做/",gene,"/7相关性图")
#vvalue<-c(pos$symbol,neg$symbol)
#for (b in 1:length(vvalue)) {
#  test <- cor.test(FF[,gene],FF[,vvalue[b]],exact=FALSE)
#  paste(paste0("n = ",length(FF[,gene])),
#        paste0("r = ",round(test$estimate,2),"(",'pearson',")"),
#        paste0("p.value= ",round(test$p.value,2)),
#        sep = ", ")
#  p<-ggplot(FF,aes(get(gene),get(vvalue[b])))+
#    geom_point(col="#984ea3")+
#    geom_smooth(method=lm, se=T,na.rm=T, fullrange=T,size=2,col="#fdc086")+
#    geom_rug(col="#7fc97f")+
#    theme_minimal()+
#    xlab(paste(gene, 'expression log2(TPM+0.001)' ,sep = ' '    ))+
#    ylab(paste(vvalue[b], 'expression log2(TPM+0.001)' ,sep = ' '    ))+
#    ## 依靠函数来生成title
#    labs(title = paste(paste0("n = ",length(FF[,gene])),
#                       paste0("r = ",round(test$estimate,2),"(",'pearson',")"),
#                      paste0("p.value= ",round(test$p.value,2)),
#                      sep = ", "))+
#   theme(plot.title = element_text(hjust = 0.5),
#          plot.margin = margin(1, 1, 1, 1, "cm"))
#  ggsave(p,filename = paste('7', gene,vvalue[b],'相关性图.pdf',sep=' '),width = 4.94,height = 4.72)
#}

FUJI<-cor_data_df %>%
  dplyr::filter(pvalue < 0.05 ) %>%
  dplyr::arrange(desc(correlation)) %>%
  top_n(300,correlation)

genename <- as.character(FUJI[,1]) #提取第一列基因名
gene_map <- biomaRt::select(org.Hs.eg.db, keys=genename, keytype="SYMBOL", columns=c("ENTREZID"))
genelist_input<-gene_map[,2]
genelist_input<-na.omit(genelist_input)
#GO分析
#BP
Go_result_BP <- enrichGO(genelist_input, 'org.Hs.eg.db', ont="BP", pvalueCutoff=1) #基因ID类型为ENSEMBL的ID形式，选择BP功能组，以P值0.05为界限
p1<-dotplot(Go_result_BP, showCategory=20)
p2<-barplot(Go_result_BP, showCategory=20)
p1<-p1 + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
p2<-p2 + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
y=Go_result_BP
yy<-as.data.frame(y)
write.csv(as.data.frame(y),paste(gene,"GO-BP.csv",sep = " "),row.names =F)#导出结果至默认路径下。
if(length(rownames(yy))!= 0) {
  ggsave(p1,filename = paste(gene,'富集分析GO_BP气泡图.pdf',sep=' '),width = 7.23,height = 8)
  ggsave(p2,filename = paste(gene,'富集分析GO_BP条形图.pdf',sep=' '),width = 7.23,height = 8)
}
#CC
Go_result_CC <- enrichGO(genelist_input, 'org.Hs.eg.db', ont="CC", pvalueCutoff=1) #基因ID类型为ENSEMBL的ID形式，选择CC功能组，以P值0.05为界限
p1<-dotplot(Go_result_CC, showCategory=20) #气泡图，显示前二十个
p2<-barplot(Go_result_CC, showCategory=20) #条形图，显示前二十个
p1<-p1 + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
p2<-p2 + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
y=Go_result_CC
yy<-as.data.frame(y)
write.csv(as.data.frame(y),paste(gene,"GO-CC.csv",sep = " "),row.names =F)#导出结果至默认路径下
if(length(rownames(yy))!= 0) {
  ggsave(p1,filename = paste(gene,'富集分析GO_CC气泡图.pdf',sep=' '),width = 7.23,height = 8)
  ggsave(p2,filename = paste(gene,'富集分析GO_CC条形图.pdf',sep=' '),width = 7.23,height = 8)
}
#MF
Go_result_MF <- enrichGO(genelist_input, 'org.Hs.eg.db',ont="MF", pvalueCutoff=1) #基因ID类型为ENSEMBL的ID形式，选择MF功能组，以P值0.05为界限
p1<-dotplot(Go_result_MF, showCategory=20) #气泡图，显示前二十个
p2<-barplot(Go_result_MF, showCategory=20) #条形图，显示前二十个
p1<-p1 + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
p2<-p2 + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
y=Go_result_MF
write.csv(as.data.frame(y),file = paste(gene,"GO-MF.csv",sep = " "),row.names =F)#导出结果至默认路径下
yy<-as.data.frame(y)
if(length(rownames(yy))!= 0){
  ggsave(p1,filename = paste(gene,'富集分析GO_MF气泡图.pdf',sep=' '),width = 7.23,height = 8)
  ggsave(p2,filename = paste(gene,'富集分析GO_MF条形图.pdf',sep=' '),width = 7.23,height = 8)
}
#KEGG分析
KEGG_result <- enrichKEGG(genelist_input, keyType = "kegg",pvalueCutoff=1,qvalueCutoff=1,pAdjustMethod = "BH", minGSSize = 5, maxGSSize = 500,organism = "hsa", use_internal_data=TRUE)  #KEGG富集分析

p1<-barplot(KEGG_result, showCategory=20)#绘制条形图
p2<-dotplot(KEGG_result, showCategory=20) #气泡图，显示前二十个
p1<-p1 + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
p2<-p2 + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
y=KEGG_result
yy<-as.data.frame(y)
write.csv(as.data.frame(y), file=paste(gene,"KEGG.csv",sep = " "),row.names =F)#导出结果至默认路径下
if(length(rownames(yy))!= 0){
  ggsave(p1,filename = paste(gene,'富集分析KEGG条形图.pdf',sep=' '),width = 7.23,height = 8)
  ggsave(p2,filename = paste(gene,'富集分析KEGG气泡图.pdf',sep=' '),width = 7.23,height = 8)
}
#GSEA分析

GSEAdata<-cor_data_df %>%
  dplyr::filter(pvalue < 0.05 ) %>%
  dplyr::arrange(desc(correlation)) %>%
  dplyr::select(symbol,correlation)
genename <- as.character(GSEAdata[,1]) #提取第一列基因名
gene_map <- biomaRt::select(org.Hs.eg.db, keys=genename, keytype="SYMBOL", columns=c("ENTREZID"))
colnames(gene_map)[1]<- 'symbol'
genelist_input<-gene_map %>%
  na.omit() %>%
  inner_join(GSEAdata,by= 'symbol') %>%
  dplyr::select(ENTREZID,correlation)
geneList = genelist_input[,2]
names(geneList) = as.character(genelist_input[,1])
geneList = sort(geneList, decreasing = TRUE)

#GSEA分析——GO
Go_gseresult <- gseGO(geneList, org.Hs.eg.db,keyType = "ENTREZID", ont="bp", nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
#GSEA分析——KEGG
KEGG_gseresult <- gseKEGG(geneList, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1,use_internal_data = T)
#GSEA分析——Reactome
Go_Reactomeresult <- gsePathway(geneList, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
write.csv (Go_gseresult,     file = paste(gene,'GSEA_GO.csv',sep=' '), row.names =F)
save(Go_gseresult,file = paste(gene,'GSEA_GO.RDATA',sep=' '))
write.csv (KEGG_gseresult,   file = paste(gene,'GSEA_KEGG.csv',sep=' '), row.names =F)
save(KEGG_gseresult,file = paste(gene,'GSEA_KEGG.RDATA',sep=' '))
write.csv (Go_Reactomeresult,file = paste(gene,'GSEA_Reactome.csv',sep=' '), row.names =F)
save(Go_Reactomeresult,file = paste(gene,'GSEA_Reactome.RDATA',sep=' '))
#波浪图
p<-ridgeplot(Go_gseresult,20) #输出前十个结果
p<-p + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
ggsave(p,filename = paste(gene,'GSEA_GO波浪图.pdf',sep=' '),width = 8,height = 8)
p<-ridgeplot(KEGG_gseresult, 20) #输出前十个结果
p<-p + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
ggsave(p,filename = paste(gene,'GSEA_KEGG波浪图.pdf',sep=' '),width = 8,height = 8)
p<-ridgeplot(Go_Reactomeresult, 20) #输出前十个结果
p<-p + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
ggsave(p,filename = paste(gene,'GSEA_Reactomeresult波浪图.pdf',sep=' '),width = 8,height = 8)
#单个图
#x=Go_gseresult@result
#x$ID<-str_replace_all(x$ID,":","_")
#number<-sum(x$p.adjust < 0.05)
#for (c in 1:number) {
#  p<-gseaplot2(Go_gseresult,c,pvalue_table = TRUE)
#  ggsave(p,filename = paste(gene,x$ID[c],'GSEA_Go.pdf',sep=' '),width = 8,height = 6)
#}
#y=KEGG_gseresult@result
#number1<-sum(y$p.adjust < 0.05)
#for (c in 1:number1) {
#  p<-gseaplot2(KEGG_gseresult,c,pvalue_table = TRUE)
#  ggsave(p,filename = paste(gene,y$ID[c],'GSEA_KEGG.pdf',sep=' '),width = 8,height = 6)
#}
#z=Go_Reactomeresult@result
#number2<-sum(z$p.adjust < 0.05)
#for (c in 1:number2) {
#  p<-gseaplot2(Go_Reactomeresult,c,pvalue_table = TRUE)
#  ggsave(p,filename = paste(gene,z$ID[c],'GSEA_Reactome.pdf',sep=' '),width = 8,height = 6)
#}















