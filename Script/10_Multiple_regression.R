rm(list = ls())



load("RS_Risk_clin.RDATA")
#load("DATA.RDATA")
clin=OSdata
# clin=read.table("clinicaldata.txt",header = T,sep = "\t")
table(clin$Risk_group)

table(clin$OS)
#加载包
library(tidyverse)
library(survminer)
library(survival)
library(survivalROC)
library(RColorBrewer)
library(forestplot)
colnames(clin)
# table(clin$N_stage)
# clin$Recurrence.metastasis[clin$Recurrence.metastasis=="NO"]=0
# clin$Recurrence.metastasis[clin$Recurrence.metastasis=="YES"]=1

colnames(clin)
clin$Risk_group[clin$Risk_group=="Low_risk"]=0
clin$Risk_group[clin$Risk_group=="High_risk"]=1



table(clin$M_stage)
clin$M_stage[clin$M_stage=="M0"]=0
clin$M_stage[clin$M_stage=="M1"]=1

table(clin$N_stage)
clin$N_stage[clin$N_stage=="N0"]=0
clin$N_stage[clin$N_stage=="N1"]=1
clin$N_stage[clin$N_stage=="N2"]=2
clin$N_stage[clin$N_stage=="N3"]=3

table(clin$T_stage)
clin$T_stage[clin$T_stage=="T1"]=1
clin$T_stage[clin$T_stage=="T2"]=2
clin$T_stage[clin$T_stage=="T3"]=1
clin$T_stage[clin$T_stage=="T4"]=2

table(clin$Gender)
clin$Gender[clin$Gender=="Female"]=0
clin$Gender[clin$Gender=="Male"]=1

table(clin$Stage)
clin$Stage[clin$Stage=="I"]=1
clin$Stage[clin$Stage=="II"]=2
clin$Stage[clin$Stage=="III"]=1
clin$Stage[clin$Stage=="IV"]=2



table(clin$Recurrence.Metastasis)
clin$Recurrence.Metastasis[clin$Recurrence.Metastasis=="NO"]=1
clin$Recurrence.Metastasis[clin$Recurrence.Metastasis=="YES"]=2




rt=clin
str(rt)
rt$Recurrence.Metastasis=as.numeric(rt$Recurrence.Metastasis)
rt$Age=as.numeric(rt$Age)
rt$M_stage=as.numeric(rt$M_stage)
rt$N_stage=as.numeric(rt$N_stage)
rt$T_stage=as.numeric(rt$T_stage)
rt$Stage=as.numeric(rt$Stage)
rt$Risk_group=as.numeric(rt$Risk_group)
rt$Gender=as.numeric(rt$Gender)
colnames(rt)
str(rt)
rt=dplyr::select(rt,ID,OS.time,OS,RS,Age,Gender,T_stage,N_stage,M_stage,Stage,Recurrence.Metastasis)
#单因素独立预后分析
uniTab=data.frame()
for(i in colnames(rt[,4:ncol(rt)])){
  cox <- coxph(Surv(OS.time, OS) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  uniTab=rbind(uniTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
write.table(uniTab,file="test.uniCox.txt",sep="\t",row.names=F,quote=F)

#多因素独立预后分析
multiCox=coxph(Surv(OS.time, OS) ~ ., data = rt[,-1])
multiCoxSum=summary(multiCox)
multiTab=data.frame()
multiTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
multiTab=cbind(id=row.names(multiTab),multiTab)
write.table(multiTab,file="test.multiCox.txt",sep="\t",row.names=F,quote=F)






############绘制森林图函数############
bioForest=function(coxFile=null,forestFile=null,forestCol=null){
  #读取输入文件
  rt <- read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  #输出图形
  pdf(file=forestFile, width = 10,height = 8)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  #绘制森林图左边的临床信息
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
  
  #绘制森林图
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, forestCol, forestCol)
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
  axis(1)
  dev.off()
}
############绘制森林图函数############

bioForest(coxFile="test.uniCox.txt",forestFile=paste0("uniForest.pdf"),forestCol="red")
bioForest(coxFile="test.multiCox.txt",forestFile=paste0("multiForest.pdf"),forestCol="red")


