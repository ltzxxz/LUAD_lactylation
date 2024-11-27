rm(list = ls())


library(RColorBrewer)
library(circlize)
library(gplots)
library(viridis)
library(oompaBase)
library(tidyverse)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor


load("RS.RDATA")
load("筛选基因表达临床.RDATA")

rs_TCGA=rs$TCGA
library(survival)
library(survminer)
cut <- surv_cutpoint(rs_TCGA,'OS.time','OS','RS')
cut
plot(cut)

## 生存分析--TCGA TCGA GSE
cat <- surv_categorize(cut)
cat1=cat
rownames(cat1)=rs_TCGA$ID

fit1 <- survfit(Surv(OS.time,OS)~RS,cat)
mytheme <- theme_survminer(font.legend = c(14,"plain", "black"),
                           font.x = c(14,"plain", "black"),
                           font.y = c(14,"plain", "black")) ## 自定义主题


ggsurvplot(fit1,cat,
           palette = 'jco',
           size=1.3,
           pval=T,
           legend.labs=c("High","Low"), 
           legend.title='Score',
           xlab="Time(years)",
           ylab='Overall survival',
           ggtheme = mytheme,
           break.time.by=2,
           conf.int=T,
           risk.table=TRUE,
           risk.table.title="",
           risk.table.height=.25)
dev.off()


rs_TCGA

rs_TCGA$Risk_group=cat1[rs_TCGA$ID,"RS"]
rs_TCGA$Risk_group=ifelse(rs_TCGA$Risk_group == "high","High_risk","Low_risk")
table(rs_TCGA$Risk_group)
rs_TCGA=rs_TCGA[,c("ID","RS","Risk_group")]

OSdata=inner_join(rs_TCGA,OSdata,by="ID")

save(OSdata,file = "RS_Risk_clin.RDATA")










###############################################################################
rt=OSdata
rownames(rt)=rt$ID
rt=rt[,-1]

colnames(rt)
#引用包
library(plyr)
library(ggplot2)
library(ggpubr)

trait="Stage"                    #临床性状

rt$OS[rt$OS==1]="Dead"
rt$OS[rt$OS==0]="Alive"
#定义临床性状的颜色
bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(unique(rt[,trait]))]

#统计高低评分组病人数目
rt1=rt[,c(trait, "Risk_group")]
rt1<-na.omit(rt1)
colnames(rt1)=c("trait", "Risk_group")
table(rt1$trait)
df=as.data.frame(table(rt1))
#计算高低评分组的百分率
df=ddply(df, .(Risk_group), transform, percent = Freq/sum(Freq) * 100)
#百分比位置
df=ddply(df, .(Risk_group), transform, pos = (cumsum(Freq) - 0.5 * Freq))
df$label=paste0(sprintf("%.0f", df$percent), "%")
df$Risk_group=factor(df$Risk_group, levels=c("Low_risk", "High_risk"))

#绘制百分率图
p=ggplot(df, aes(x = factor(Risk_group), y = percent, fill = trait)) +
  geom_bar(position = position_stack(), stat = "identity", width = .7) +
  scale_fill_manual(values=bioCol)+
  xlab("score")+ ylab("Percent weight")+  guides(fill=guide_legend(title=trait))+
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  #coord_flip()+
  theme_bw()
pdf(file="Stage_barplot.pdf", width=4, height=5)
print(p)
dev.off()

#设置比较组
rt2=rt[,c(trait, "RS")]
rt2<-na.omit(rt2)
colnames(rt2)=c("trait", "score")
type=levels(factor(rt2[,"trait"]))
rt2$trait=factor(rt2$trait,levels =type )
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
#绘制箱线图
boxplot=ggboxplot(rt2, x="trait", y="score", fill="trait",
                  xlab=trait,
                  ylab="score",
                  legend.title=trait,
                  palette=bioCol
)+ 
  stat_compare_means(comparisons=my_comparisons)
pdf(file="Stage_boxplot.pdf",width=5,height=4.5)
print(boxplot)
dev.off()

colnames(rt)

