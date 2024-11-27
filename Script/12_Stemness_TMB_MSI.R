rm(list = ls())
library(tidyverse)
library(ComplexHeatmap)
library(cowplot)
library(ggpubr)
load("RS_Risk_clin.RDATA")
clin=OSdata
group=clin[,c("ID","RS", "Risk_group")]
mRNAsi=read.table("泛癌干性mRNAsi.txt",header = T,sep = "\t")
mDNAsi=read.table("泛癌干性mDNAsi.txt",header = T,sep = "\t")
TMB=read.table("泛癌TMB.txt",header = T,sep = "\t")
MSI=read.table("泛癌MSI.txt",header = T,sep = "\t")
colnames(group)[1]="ID"
colnames(mRNAsi)[1]="ID"
colnames(mDNAsi)[1]="ID"
colnames(TMB)[1]="ID"
colnames(MSI)[1]="ID"
colnames(clin)[1]="ID"
aaa=intersect(clin$ID,group$ID)
clin=dplyr::filter(clin,ID %in% aaa)

groupdata=clin
# groupdata=dplyr::filter(groupdata,cohort == "TCGA")
# groupdata$ID=str_replace_all(groupdata$ID,"TCGA_702_","")
# groupdata$ID=str_replace_all(groupdata$ID,"_",".")
#mRNAsi
mRNAsi$ID=substr(mRNAsi$ID,1,15)
mRNAsi$ID=str_replace_all(mRNAsi$ID,"-",".")
bbb=intersect(groupdata$ID,mRNAsi$ID)
group1=dplyr::filter(groupdata,ID %in% bbb)
mRNAsi=dplyr::filter(mRNAsi,ID %in% bbb)
mRNAsidata=inner_join(group1,mRNAsi,by="ID")


input=mRNAsidata
input=dplyr::select(input,ID,Risk_group,OS,mRNAsi)
input <- as.matrix(input)
input[is.na(input)] <- "Unknown"
input <- as.data.frame(input)
input$mRNAsi <- as.numeric(input$mRNAsi)


input$Risk_group <- factor(input$Risk_group,levels = c("High_risk","Low_risk"))
input$OS[input$OS==1]="Dead"
input$OS[input$OS==0]="Alive"
input$OS <- factor(input$OS,levels = c("Dead","Alive"))
input=arrange(input,mRNAsi)
input$index <- 1:nrow(input)

# 画图
darkred   <- "#F2042C"
blue      <- "#0A0745"
lightgrey <- "#dcddde"
# 自定义主题
My_Theme1 <- theme_bw() + 
  theme(panel.grid =element_blank()) + 
  theme(panel.border = element_rect(size = 1)) + 
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size = 12,angle = 90),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
  )+
  theme(plot.margin = margin(0,0.1,0,0.1, "cm"))

My_Theme2 = theme_minimal()+
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size = 12,angle = 0),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  theme (panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.margin = margin(0,0.1,0,0.1, "cm"))

p1 <- ggplot(input, aes(index, mRNAsi))+
  geom_bar(stat="identity",col = blue)+
  My_Theme1 +
  labs(y = "mRNAsi") +
  scale_x_continuous(expand = c(0,0)) #不留空

p2 <- ggplot(input,aes(index,1))+
  geom_tile(aes(fill = OS))+
  My_Theme2+
  labs(y = "OS")+
  scale_fill_manual(values = c("#74C065","#B096C2","#FAA96C","#EFF882","#6D6466")) +
  scale_x_continuous(expand = c(0,0)) #不留空

p3 <- ggplot(input,aes(index,1))+
  geom_tile(aes(fill = Risk_group))+
  My_Theme2+
  labs(y = "Risk_group")+
  scale_fill_manual(values = c("#84CBA8","#F2FA9C","#B1A5C8","#FA5E5C")) +
  scale_x_continuous(expand = c(0,0)) #不留空


legend_a <- get_legend(p2+theme(legend.position = "bottom"))
legend_b <- get_legend(p3+theme(legend.position = "bottom"))

p <- plot_grid(p1,p2,p3,
               align = "v",axis = "l",
               ncol = 1, rel_heights = c(4,1,1,1,1),
               legend_a,legend_b)
p

save_plot("mRNAsi.pdf", p,base_width = 10,base_height = 8)


data=mRNAsidata
#设置比较组
data$Risk_group=factor(data$Risk_group, levels=c("Low_risk", "High_risk"))
my_comparisons <- list( c("Low_risk", "High_risk") )

data$mRNAsi=as.numeric(data$mRNAsi)
#绘制boxplot


p <- ggboxplot(data, x = "Risk_group", y = "mRNAsi",
               fill = "Risk_group",legend=F,palette =c("#E7B800", "#00AFBB"),bxp.errorbar=T)+
  theme(legend.position='none')+
  stat_compare_means(comparisons = my_comparisons,method = 't.test',aes(label = ..p.signif..))
ggsave(p,filename = paste("mRNAsi",'差异表达.pdf',sep=' '),width = 5.6,height = 4.22)

#相关性分析、
inputtemp=mRNAsidata
test <- cor.test(inputtemp[,"RS"],inputtemp[,"mRNAsi"],exact=FALSE)
paste(paste0("n = ",length(inputtemp[,"RS"])),
      paste0("r = ",round(test$estimate,2),"(",'pearson',")"),
      paste0("p.value= ",round(test$p.value,4)),
      sep = ", ")
p<-ggplot(inputtemp,aes(RS,mRNAsi))+
  geom_point(col="#984ea3")+
  geom_smooth(method=lm, se=T,na.rm=T, fullrange=T,size=2,col="#fdc086")+
  geom_rug(col="#7fc97f")+
  theme_minimal()+
  xlab(paste("score" ,sep = ' '    ))+
  ylab(paste("mRNAsi",sep = ' '    ))+
  ## 依靠函数来生成title
  labs(title = paste(paste0("n = ",length(inputtemp[,"RS"])),
                     paste0("r = ",round(test$estimate,2),"(",'pearson',")"),
                     paste0("p.value= ",round(test$p.value,4)),
                     sep = ", "))+
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin = margin(1, 1, 1, 1, "cm"))
ggsave(p,filename = paste("mRNAsi","score",'相关性点图.pdf',sep=' '),width = 4.94,height = 4.72)
dev.off()

#####################################################################################################
#mDNAsi
mDNAsi$ID=substr(mDNAsi$ID,1,15)
mDNAsi$ID=str_replace_all(mDNAsi$ID,"-",".")
bbb=intersect(groupdata$ID,mDNAsi$ID)
group1=dplyr::filter(groupdata,ID %in% bbb)
mDNAsi=dplyr::filter(mDNAsi,ID %in% bbb)
mDNAsidata=inner_join(group1,mDNAsi,by="ID")


input=mDNAsidata
input=dplyr::select(input,ID,Risk_group,OS,mDNAsi)
input <- as.matrix(input)
input[is.na(input)] <- "Unknown"
input <- as.data.frame(input)
input$mDNAsi <- as.numeric(input$mDNAsi)


input$Risk_group <- factor(input$Risk_group,levels = c("High_risk","Low_risk"))
input$OS[input$OS==1]="Dead"
input$OS[input$OS==0]="Alive"
input$OS <- factor(input$OS,levels = c("Dead","Alive"))
input=arrange(input,mDNAsi)
input$index <- 1:nrow(input)

# 画图
darkred   <- "#F2042C"
blue      <- "#0A0745"
lightgrey <- "#dcddde"
# 自定义主题
My_Theme1 <- theme_bw() + 
  theme(panel.grid =element_blank()) + 
  theme(panel.border = element_rect(size = 1)) + 
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size = 12,angle = 90),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
  )+
  theme(plot.margin = margin(0,0.1,0,0.1, "cm"))

My_Theme2 = theme_minimal()+
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size = 12,angle = 0),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  theme (panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.margin = margin(0,0.1,0,0.1, "cm"))

p1 <- ggplot(input, aes(index, mDNAsi))+
  geom_bar(stat="identity",col = blue)+
  My_Theme1 +
  labs(y = "mDNAsi") +
  scale_x_continuous(expand = c(0,0)) #不留空

p2 <- ggplot(input,aes(index,1))+
  geom_tile(aes(fill = OS))+
  My_Theme2+
  labs(y = "OS")+
  scale_fill_manual(values = c("#74C065","#B096C2","#FAA96C","#EFF882","#6D6466")) +
  scale_x_continuous(expand = c(0,0)) #不留空

p3 <- ggplot(input,aes(index,1))+
  geom_tile(aes(fill = Risk_group))+
  My_Theme2+
  labs(y = "group")+
  scale_fill_manual(values = c("#84CBA8","#F2FA9C","#B1A5C8","#FA5E5C")) +
  scale_x_continuous(expand = c(0,0)) #不留空


legend_a <- get_legend(p2+theme(legend.position = "bottom"))
legend_b <- get_legend(p3+theme(legend.position = "bottom"))

p <- plot_grid(p1,p2,p3,
               align = "v",axis = "l",
               ncol = 1, rel_heights = c(4,1,1,1,1),
               legend_a,legend_b)
p

save_plot("mDNAsi.pdf", p,base_width = 10,base_height = 8)


data=mDNAsidata
#设置比较组
data$Risk_group=factor(data$Risk_group, levels=c("Low_risk", "High_risk"))
my_comparisons <- list( c("Low_risk", "High_risk") )

data$mDNAsi=as.numeric(data$mDNAsi)
#绘制boxplot


p <- ggboxplot(data, x = "Risk_group", y = "mDNAsi",
               fill = "Risk_group",legend=F,palette =c("#E7B800", "#00AFBB"),bxp.errorbar=T)+
  theme(legend.position='none')+
  stat_compare_means(comparisons = my_comparisons,method = 't.test',aes(label = ..p.signif..))
ggsave(p,filename = paste("mDNAsi",'差异表达.pdf',sep=' '),width = 5.6,height = 4.22)

#相关性分析、
inputtemp=mDNAsidata
test <- cor.test(inputtemp[,"RS"],inputtemp[,"mDNAsi"],exact=FALSE)
paste(paste0("n = ",length(inputtemp[,"RS"])),
      paste0("r = ",round(test$estimate,2),"(",'pearson',")"),
      paste0("p.value= ",round(test$p.value,4)),
      sep = ", ")
p<-ggplot(inputtemp,aes(RS,mDNAsi))+
  geom_point(col="#984ea3")+
  geom_smooth(method=lm, se=T,na.rm=T, fullrange=T,size=2,col="#fdc086")+
  geom_rug(col="#7fc97f")+
  theme_minimal()+
  xlab(paste("score" ,sep = ' '    ))+
  ylab(paste("mDNAsi",sep = ' '    ))+
  ## 依靠函数来生成title
  labs(title = paste(paste0("n = ",length(inputtemp[,"RS"])),
                     paste0("r = ",round(test$estimate,2),"(",'pearson',")"),
                     paste0("p.value= ",round(test$p.value,4)),
                     sep = ", "))+
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin = margin(1, 1, 1, 1, "cm"))
ggsave(p,filename = paste("mDNAsi","score",'相关性点图.pdf',sep=' '),width = 4.94,height = 4.72)


#####################################################################################################
#TMB
TMB$ID=substr(TMB$ID,1,15)
TMB$ID=str_replace_all(TMB$ID,"-",".")
groupdata1=groupdata
groupdata1$ID=substr(groupdata1$ID,1,15)
bbb=intersect(groupdata1$ID,TMB$ID)
group1=dplyr::filter(groupdata1,ID %in% bbb)
TMB=dplyr::filter(TMB,ID %in% bbb)
TMBdata=inner_join(group1,TMB,by="ID")


input=TMBdata
input=dplyr::select(input,ID,Risk_group,OS,TMB)
input <- as.matrix(input)
input[is.na(input)] <- "Unknown"
input <- as.data.frame(input)
input$TMB <- as.numeric(input$TMB)


input$Risk_group <- factor(input$Risk_group,levels = c("High_risk","Low_risk"))
input$OS[input$OS==1]="Dead"
input$OS[input$OS==0]="Alive"
input$OS <- factor(input$OS,levels = c("Dead","Alive"))
input=arrange(input,TMB)
input$index <- 1:nrow(input)

# 画图
darkred   <- "#F2042C"
blue      <- "#0A0745"
lightgrey <- "#dcddde"
# 自定义主题
My_Theme1 <- theme_bw() + 
  theme(panel.grid =element_blank()) + 
  theme(panel.border = element_rect(size = 1)) + 
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size = 12,angle = 90),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
  )+
  theme(plot.margin = margin(0,0.1,0,0.1, "cm"))

My_Theme2 = theme_minimal()+
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size = 12,angle = 0),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  theme (panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.margin = margin(0,0.1,0,0.1, "cm"))

p1 <- ggplot(input, aes(index, TMB))+
  geom_bar(stat="identity",col = blue)+
  My_Theme1 +
  labs(y = "TMB") +
  scale_x_continuous(expand = c(0,0)) #不留空

p2 <- ggplot(input,aes(index,1))+
  geom_tile(aes(fill = OS))+
  My_Theme2+
  labs(y = "OS")+
  scale_fill_manual(values = c("#74C065","#B096C2","#FAA96C","#EFF882","#6D6466")) +
  scale_x_continuous(expand = c(0,0)) #不留空

p3 <- ggplot(input,aes(index,1))+
  geom_tile(aes(fill = Risk_group))+
  My_Theme2+
  labs(y = "group")+
  scale_fill_manual(values = c("#84CBA8","#F2FA9C","#B1A5C8","#FA5E5C")) +
  scale_x_continuous(expand = c(0,0)) #不留空


legend_a <- get_legend(p2+theme(legend.position = "bottom"))
legend_b <- get_legend(p3+theme(legend.position = "bottom"))

p <- plot_grid(p1,p2,p3,
               align = "v",axis = "l",
               ncol = 1, rel_heights = c(4,1,1,1,1),
               legend_a,legend_b)
p

save_plot("TMB.pdf", p,base_width = 10,base_height = 8)


data=TMBdata
#设置比较组
data$Risk_group=factor(data$Risk_group, levels=c("Low_risk", "High_risk"))
my_comparisons <- list( c("Low_risk", "High_risk") )

data$TMB=as.numeric(data$TMB)
#绘制boxplot


p <- ggboxplot(data, x = "Risk_group", y = "TMB",
               fill = "Risk_group",legend=F,palette =c("#E7B800", "#00AFBB"),bxp.errorbar=T)+
  theme(legend.position='none')+
  stat_compare_means(comparisons = my_comparisons,method = 't.test',aes(label = ..p.signif..))
ggsave(p,filename = paste("TMB",'差异表达.pdf',sep=' '),width = 5.6,height = 4.22)

#相关性分析、
inputtemp=TMBdata
test <- cor.test(inputtemp[,"RS"],inputtemp[,"TMB"],exact=FALSE)
paste(paste0("n = ",length(inputtemp[,"RS"])),
      paste0("r = ",round(test$estimate,2),"(",'pearson',")"),
      paste0("p.value= ",round(test$p.value,4)),
      sep = ", ")
p<-ggplot(inputtemp,aes(RS,TMB))+
  geom_point(col="#984ea3")+
  geom_smooth(method=lm, se=T,na.rm=T, fullrange=T,size=2,col="#fdc086")+
  geom_rug(col="#7fc97f")+
  theme_minimal()+
  xlab(paste("score" ,sep = ' '    ))+
  ylab(paste("TMB",sep = ' '    ))+
  ## 依靠函数来生成title
  labs(title = paste(paste0("n = ",length(inputtemp[,"RS"])),
                     paste0("r = ",round(test$estimate,2),"(",'pearson',")"),
                     paste0("p.value= ",round(test$p.value,4)),
                     sep = ", "))+
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin = margin(1, 1, 1, 1, "cm"))
ggsave(p,filename = paste("TMB","score",'相关性点图.pdf',sep=' '),width = 4.94,height = 4.72)


#####################################################################################################
#MSI
MSI$ID=substr(MSI$ID,1,15)
MSI$ID=str_replace_all(MSI$ID,"-",".")
groupdata1=groupdata
groupdata1$ID=substr(groupdata1$ID,1,16)
bbb=intersect(groupdata1$ID,MSI$ID)
group1=dplyr::filter(groupdata1,ID %in% bbb)
MSI=dplyr::filter(MSI,ID %in% bbb)
MSIdata=inner_join(group1,MSI,by="ID")


input=MSIdata
input=dplyr::select(input,ID,Risk_group,OS,MSI)
input <- as.matrix(input)
input[is.na(input)] <- "Unknown"
input <- as.data.frame(input)
input$MSI <- as.numeric(input$MSI)


input$Risk_group <- factor(input$Risk_group,levels = c("High_risk","Low_risk"))
input$OS[input$OS==1]="Dead"
input$OS[input$OS==0]="Alive"
input$OS <- factor(input$OS,levels = c("Dead","Alive"))
input=arrange(input,MSI)
input$index <- 1:nrow(input)

# 画图
darkred   <- "#F2042C"
blue      <- "#0A0745"
lightgrey <- "#dcddde"
# 自定义主题
My_Theme1 <- theme_bw() + 
  theme(panel.grid =element_blank()) + 
  theme(panel.border = element_rect(size = 1)) + 
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size = 12,angle = 90),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
  )+
  theme(plot.margin = margin(0,0.1,0,0.1, "cm"))

My_Theme2 = theme_minimal()+
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size = 12,angle = 0),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  theme (panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.margin = margin(0,0.1,0,0.1, "cm"))

p1 <- ggplot(input, aes(index, MSI))+
  geom_bar(stat="identity",col = blue)+
  My_Theme1 +
  labs(y = "MSI") +
  scale_x_continuous(expand = c(0,0)) #不留空

p2 <- ggplot(input,aes(index,1))+
  geom_tile(aes(fill = OS))+
  My_Theme2+
  labs(y = "OS")+
  scale_fill_manual(values = c("#74C065","#B096C2","#FAA96C","#EFF882","#6D6466")) +
  scale_x_continuous(expand = c(0,0)) #不留空

p3 <- ggplot(input,aes(index,1))+
  geom_tile(aes(fill = Risk_group))+
  My_Theme2+
  labs(y = "group")+
  scale_fill_manual(values = c("#84CBA8","#F2FA9C","#B1A5C8","#FA5E5C")) +
  scale_x_continuous(expand = c(0,0)) #不留空


legend_a <- get_legend(p2+theme(legend.position = "bottom"))
legend_b <- get_legend(p3+theme(legend.position = "bottom"))

p <- plot_grid(p1,p2,p3,
               align = "v",axis = "l",
               ncol = 1, rel_heights = c(4,1,1,1,1),
               legend_a,legend_b)
p

save_plot("MSI.pdf", p,base_width = 10,base_height = 8)


data=MSIdata
#设置比较组
data$Risk_group=factor(data$Risk_group, levels=c("Low_risk", "High_risk"))
my_comparisons <- list( c("Low_risk", "High_risk") )

data$MSI=as.numeric(data$MSI)
#绘制boxplot


p <- ggboxplot(data, x = "Risk_group", y = "MSI",
               fill = "Risk_group",legend=F,palette =c("#E7B800", "#00AFBB"),bxp.errorbar=T)+
  theme(legend.position='none')+
  stat_compare_means(comparisons = my_comparisons,method = 't.test',aes(label = ..p.signif..))
ggsave(p,filename = paste("MSI",'差异表达.pdf',sep=' '),width = 5.6,height = 4.22)

#相关性分析、
inputtemp=MSIdata
test <- cor.test(inputtemp[,"RS"],inputtemp[,"MSI"],exact=FALSE)
paste(paste0("n = ",length(inputtemp[,"RS"])),
      paste0("r = ",round(test$estimate,2),"(",'pearson',")"),
      paste0("p.value= ",round(test$p.value,4)),
      sep = ", ")
p<-ggplot(inputtemp,aes(RS,MSI))+
  geom_point(col="#984ea3")+
  geom_smooth(method=lm, se=T,na.rm=T, fullrange=T,size=2,col="#fdc086")+
  geom_rug(col="#7fc97f")+
  theme_minimal()+
  xlab(paste("score" ,sep = ' '    ))+
  ylab(paste("MSI",sep = ' '    ))+
  ## 依靠函数来生成title
  labs(title = paste(paste0("n = ",length(inputtemp[,"RS"])),
                     paste0("r = ",round(test$estimate,2),"(",'pearson',")"),
                     paste0("p.value= ",round(test$p.value,4)),
                     sep = ", "))+
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin = margin(1, 1, 1, 1, "cm"))
ggsave(p,filename = paste("MSI","score",'相关性点图.pdf',sep=' '),width = 4.94,height = 4.72)


#####################################################################################################










