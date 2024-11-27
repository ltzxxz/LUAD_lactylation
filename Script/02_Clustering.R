rm(list = ls())
library(tidyverse)
library(data.table)
library(maftools)
library("MOVICS")

load("cleandata.RDATA")


mo.data   <- brca.tcga[1:4]
genename=read.table("Lactylation.txt",header = T)
genename=intersect(genename$ID,rownames(mo.data$mRNA.expr))
mo.data$mRNA.expr=mo.data$mRNA.expr[genename,]
brca.tcga$mRNA.expr=mo.data$mRNA.expr
# scenario 4: considering we are dealing with data and use cox to select elite
surv.info <- brca.tcga$clin.info
tmp       <- brca.tcga$mRNA.expr # get expression data 
elite.tmp <- getElites(dat       = tmp,
                       method    = "cox",
                       surv.info = surv.info, # must provide survival information with 'futime' and 'fustat'
                       p.cutoff  = 0.05,
                       elite.num = 100) # this time elite.num argument will be disabled because cox method refers to p.cutoff to select elites

dim(elite.tmp$elite.dat) 

table(elite.tmp$unicox$pvalue < 0.01) #  genes have nominal pvalue < 0.05 in univariate Cox regression
coxdata=elite.tmp$unicox
coxdata$pvalue=as.numeric(coxdata$pvalue)
coxdata=dplyr::filter(coxdata,pvalue < 0.01)
geneid=coxdata$gene
mo.data$mRNA.expr=mo.data$mRNA.expr[geneid,]

tmp       <- brca.tcga$lcRNA.expr # get expression data 
elite.tmp <- getElites(dat       = tmp,
                       method    = "cox",
                       surv.info = surv.info, # must provide survival information with 'futime' and 'fustat'
                       p.cutoff  = 0.05,
                       elite.num = 100) # this time elite.num argument will be disabled because cox method refers to p.cutoff to select elites

dim(elite.tmp$elite.dat) 

table(elite.tmp$unicox$pvalue < 0.01) #  genes have nominal pvalue < 0.05 in univariate Cox regression
coxdata=elite.tmp$unicox
coxdata$pvalue=as.numeric(coxdata$pvalue)
coxdata=dplyr::filter(coxdata,pvalue < 0.01)
geneid=coxdata$gene
mo.data$lcRNA.expr=mo.data$lcRNA.expr[geneid,]

tmp       <- brca.tcga$meth.beta # get expression data 
tmp=na.omit(tmp)
str(tmp)

elite.tmp <- getElites(dat       = tmp,
                       method    = "cox",
                       surv.info = surv.info, # must provide survival information with 'futime' and 'fustat'
                       p.cutoff  = 0.05,
                       elite.num = 100) # this time elite.num argument will be disabled because cox method refers to p.cutoff to select elites

dim(elite.tmp$elite.dat) 

table(elite.tmp$unicox$pvalue < 0.01) #  genes have nominal pvalue < 0.05 in univariate Cox regression
coxdata=elite.tmp$unicox
coxdata$pvalue=as.numeric(coxdata$pvalue)
coxdata=dplyr::filter(coxdata,pvalue < 0.01)
geneid=coxdata$gene
rownames(mo.data$meth.beta)
mo.data$meth.beta=mo.data$meth.beta[geneid,]



tmp       <- brca.tcga$mut.status # get expression data 
elite.tmp <- getElites(dat       = tmp,
                       method    = "cox",
                       surv.info = surv.info, # must provide survival information with 'futime' and 'fustat'
                       p.cutoff  = 0.05,
                       elite.num = 100) # this time elite.num argument will be disabled because cox method refers to p.cutoff to select elites

dim(elite.tmp$elite.dat) 

table(elite.tmp$unicox$pvalue < 0.01) #  genes have nominal pvalue < 0.05 in univariate Cox regression
coxdata=elite.tmp$unicox
coxdata$pvalue=as.numeric(coxdata$pvalue)
coxdata=dplyr::filter(coxdata,pvalue < 0.01)
geneid=coxdata$gene
mo.data$mut.status=mo.data$mut.status[geneid,]







# identify optimal clustering number (may take a while)
optk.brca <- getClustNum(data        = mo.data,
                         is.binary   = c(F,F,F,T), # note: the 4th data is somatic mutation which is a binary matrix
                         try.N.clust = 2:8, # try cluster number from 2 to 8
                         fig.name    = "CLUSTER NUMBER OF TCGA-LUAD")




# perform iClusterBayes (may take a while)


iClusterBayes.res <- getMOIC(data        = mo.data,
                             N.clust     = 2,
                             methodslist = "iClusterBayes", # specify only ONE algorithm here
                             type        = c("gaussian","gaussian","gaussian","binomial"), # data type corresponding to the list
                             n.burnin    = 1800,
                             n.draw      = 1200,
                             prior.gamma = c(0.5, 0.5, 0.5, 0.5),
                             sdev        = 0.05,
                             thin        = 3)


# perform multi-omics integrative clustering with the rest of 9 algorithms
moic.res.list <- getMOIC(data        = mo.data,
                         methodslist = list("SNF", "PINSPlus", "NEMO", "COCA", "LRAcluster", "ConsensusClustering", "IntNMF", "CIMLR", "MoCluster"),
                         N.clust     = 2,
                         type        = c("gaussian", "gaussian", "gaussian", "binomial"))
moic.res.list <- append(moic.res.list, 
                        list("iClusterBayes" = iClusterBayes.res))
# save moic.res.list to local path
save(moic.res.list, file = "moic.res.list.rda")





cmoic.brca <- getConsensusMOIC(moic.res.list = moic.res.list,
                               fig.name      = "CONSENSUS HEATMAP",
                               distance      = "euclidean",
                               linkage       = "average")

cluster=cmoic.brca$clust.res
colnames(cluster)=c("ID",	"cluster")
cluster$cluster[cluster$cluster == 1]="CS1"
cluster$cluster[cluster$cluster == 2]="CS2"
write.table(cluster,file = "cluster.txt",sep = "\t",quote = F,row.names = F)


getSilhouette(sil      = cmoic.brca$sil, # a sil object returned by getConsensusMOIC()
              fig.path = getwd(),
              fig.name = "SILHOUETTE",
              height   = 5.5,
              width    = 5)

# convert beta value to M value for stronger signal
indata <- mo.data
# indata$meth.beta <- log2(indata$meth.beta / (1 - indata$meth.beta))

# data normalization for heatmap
plotdata <- getStdiz(data       = indata,
                     halfwidth  = c(2,2,2,NA), # no truncation for mutation
                     centerFlag = c(T,T,T,F), # no center for mutation
                     scaleFlag  = c(T,T,T,F)) # no scale for mutation


feat   <- iClusterBayes.res$feat.res
feat1  <- feat[which(feat$dataset == "mRNA.expr"),][1:10,"feature"] 
feat2  <- feat[which(feat$dataset == "lcRNA.expr"),][1:10,"feature"]
feat3  <- feat[which(feat$dataset == "meth.beta"),][1:10,"feature"]
feat4  <- feat[which(feat$dataset == "mut.status"),][1:10,"feature"]
annRow <- list(feat1, feat2, feat3, feat4)



# set color for each omics data
# if no color list specified all subheatmaps will be unified to green and red color pattern
mRNA.col   <- c("#00FF00", "#008000", "#000000", "#800000", "#FF0000")
lcRNA.col <- c("#6699CC", "white"  , "#FF3C38")
meth.col   <- c("#0074FE", "#96EBF9", "#FEE900", "#F00003")
mut.col    <- c("grey90" , "black")
col.list   <- list(mRNA.col, lcRNA.col, meth.col, mut.col)

# comprehensive heatmap (may take a while)
getMoHeatmap(data          = plotdata,
             row.title     = c("mRNA","lcRNA","Methylation","Mutation"),
             is.binary     = c(F,F,F,T), # the 4th data is mutation which is binary
             legend.name   = c("mRNA.FPKM","lcRNA.FPKM","M value","Mutated"),
             clust.res     = iClusterBayes.res$clust.res, # cluster results
             clust.dend    = NULL, # no dendrogram
             show.rownames = c(F,F,F,F), # specify for each omics data
             show.colnames = FALSE, # show no sample names
             annRow        = annRow, # mark selected features
             color         = col.list,
             annCol        = NULL, # no annotation for samples
             annColors     = NULL, # no annotation color
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "COMPREHENSIVE HEATMAP OF ICLUSTERBAYES")



# comprehensive heatmap (may take a while)
getMoHeatmap(data          = plotdata,
             row.title     = c("mRNA","lcRNA","Methylation","Mutation"),
             is.binary     = c(F,F,F,T), # the 4th data is mutation which is binary
             legend.name   = c("mRNA.FPKM","lcRNA.FPKM","M value","Mutated"),
             clust.res     = moic.res.list$COCA$clust.res, # cluster results
             clust.dend    = moic.res.list$COCA$clust.dend, # show dendrogram for samples
             color         = col.list,
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "COMPREHENSIVE HEATMAP OF COCA")
#####生存分析
# survival comparison
surv.brca <- compSurv(moic.res         = cmoic.brca,
                      surv.info        = surv.info,
                      convt.time       = "y", # convert day unit to month
                      surv.median.line = "h", # draw horizontal line at median survival
                      xyrs.est         = c(5,10), # estimate 5 and 10-year survival
                      fig.name         = "KAPLAN-MEIER CURVE OF CONSENSUSMOIC")


pdf(file = "KAPLAN-MEIER CURVE.pdf",width = 8,height = 6)
compSurv(moic.res         = cmoic.brca,
         surv.info        = surv.info,
         convt.time       = "y", # convert day unit to month
         surv.median.line = "h", # draw horizontal line at median survival
         xyrs.est         = c(5,10), # estimate 5 and 10-year survival
         fig.name         = "KAPLAN-MEIER CURVE OF CONSENSUSMOIC")
dev.off()



# mutational frequency comparison
aaa=brca.tcga$mut.status
aaa$sums <- rowSums(aaa)
aaa <- aaa[order(aaa$sums, decreasing = TRUE), ]
genes=rownames(aaa)[1:20]
mut.brca <- compMut(moic.res     = cmoic.brca,
                    mut.matrix   = brca.tcga$mut.status[genes,], # binary somatic mutation matrix
                    doWord       = TRUE, # generate table in .docx format
                    doPlot       = TRUE, # draw OncoPrint
                    freq.cutoff  = 0.05, # keep those genes that mutated in at least 5% of samples
                    p.adj.cutoff = 0.05, # keep those genes with adjusted p value < 0.05 to draw OncoPrint
                    innerclust   = TRUE, # perform clustering within each subtype
                    #annCol       = annCol, # same annotation for heatmap
                    #annColors    = annColors, # same annotation color for heatmap
                    width        = 6, 
                    height       = 4,
                    fig.name     = "ONCOPRINT FOR SIGNIFICANT MUTATIONS",
                    tab.name     = "INDEPENDENT TEST BETWEEN SUBTYPE AND MUTATION")
print(mut.brca)

write.table(mut.brca,file = "Comparison of mutational frequency.txt",sep = "\t",quote = F,row.names = F)

###差异分析
mRNA.expr=mo.data$mRNA.expr
runDEA(dea.method = "limma",
       expr       = mRNA.expr, # normalized expression data
       moic.res   = cmoic.brca,
       prefix     = "TCGA-LUAD")

# choose limma result to identify subtype-specific up-regulated biomarkers
marker.up <- runMarker(moic.res      = cmoic.brca,
                       dea.method    = "limma", # name of DEA method
                       prefix        = "TCGA-LUAD", # MUST be the same of argument in runDEA()
                       dat.path      = getwd(), # path of DEA files
                       res.path      = getwd(), # path to save marker files
                       p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                       p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                       dirct         = "up", # direction of dysregulation in expression
                       n.marker      = 100, # number of biomarkers for each subtype
                       doplot        = TRUE, # generate diagonal heatmap
                       norm.expr     = mRNA.expr, # use normalized expression as heatmap input
                       #annCol        = annCol, # sample annotation in heatmap
                       #annColors     = annColors, # colors for sample annotation
                       show_rownames = FALSE, # show no rownames (biomarker name)
                       fig.name      = "UPREGULATED BIOMARKER HEATMAP")
#> --all samples matched.
#> --log2 transformation done for expression data.
# check the upregulated biomarkers
head(marker.up$templates)

MSIGDB.FILE <- system.file("extdata", "c5.bp.v7.1.symbols.xls", package = "MOVICS", mustWork = TRUE)


#run GSEA to identify up-regulated GO pathways using results from edgeR
gsea.up <- runGSEA(moic.res     = cmoic.brca,
                   dea.method   = "limma", # name of DEA method
                   prefix       = "TCGA-LUAD", # MUST be the same of argument in runDEA()
                   dat.path     = getwd(), # path of DEA files
                   res.path     = getwd(), # path to save GSEA files
                   msigdb.path  = MSIGDB.FILE, # MUST be the ABSOLUTE path of msigdb file
                   norm.expr    = mRNA.expr, # use normalized expression to calculate enrichment score
                   dirct        = "up", # direction of dysregulation in pathway
                   p.cutoff     = 0.05, # p cutoff to identify significant pathways
                   p.adj.cutoff = 0.25, # padj cutoff to identify significant pathways
                   gsva.method  = "gsva", # method to calculate single sample enrichment score
                   norm.method  = "mean", # normalization method to calculate subtype-specific enrichment score
                   fig.name     = "UPREGULATED PATHWAY HEATMAP")



# MUST locate ABSOLUTE path of gene set file
GSET.FILE <- 
  system.file("extdata", "gene sets of interest.gmt", package = "MOVICS", mustWork = TRUE)

# run GSVA to estimate single sample enrichment score based on given gene set of interest
gsva.res <- 
  runGSVA(moic.res      = cmoic.brca,
          norm.expr     = mRNA.expr,
          gset.gmt.path = GSET.FILE, # ABSOLUTE path of gene set file
          gsva.method   = "gsva", # method to calculate single sample enrichment score
          annCol        = annCol,
          #annColors     = annColors,
          fig.path      = getwd(),
          fig.name      = "GENE SETS OF INTEREST HEATMAP",
          height        = 5,
          width         = 8)

######4) run nearest template prediction in external cohort
load("GSE13213.RDATA")
# run NTP in Yau cohort by using up-regulated biomarkers

yau.ntp.pred <- runNTP(expr       = GSE13213$mRNA.expr,
                       templates  = marker.up$templates, # the template has been already prepared in runMarker()
                       scaleFlag  = TRUE, # scale input data (by default)
                       centerFlag = TRUE, # center input data (by default)
                       doPlot     = TRUE, # to generate heatmap
                       fig.name   = "GSE13213") 

# compare survival outcome in Yau cohort
surv.yau <- compSurv(moic.res         = yau.ntp.pred,
                     surv.info        = GSE13213$clin.info,
                     convt.time       = "y", # switch to month
                     surv.median.line = "hv", # switch to both
                     fig.name         = "KAPLAN-MEIER CURVE OF GSE13213") 
dev.off()
pdf(file = "KAPLAN-MEIER CURVE OF GSE13213.pdf",width = 8,height = 6)
compSurv(moic.res         = yau.ntp.pred,
         surv.info        = GSE13213$clin.info,
         convt.time       = "y", # switch to month
         surv.median.line = "hv", # switch to both
         fig.name         = "KAPLAN-MEIER CURVE OF GSE13213") 
dev.off()


######4) run nearest template prediction in external cohort
load("GSE31210.RDATA")
# run NTP in Yau cohort by using up-regulated biomarkers
yau.ntp.pred <- runNTP(expr       = GSE31210$mRNA.expr,
                       templates  = marker.up$templates, # the template has been already prepared in runMarker()
                       scaleFlag  = TRUE, # scale input data (by default)
                       centerFlag = TRUE, # center input data (by default)
                       doPlot     = TRUE, # to generate heatmap
                       fig.name   = "GSE31210") 

# compare survival outcome in Yau cohort
surv.yau <- compSurv(moic.res         = yau.ntp.pred,
                     surv.info        = GSE31210$clin.info,
                     convt.time       = "y", # switch to month
                     surv.median.line = "hv", # switch to both
                     fig.name         = "KAPLAN-MEIER CURVE OF GSE31210") 
dev.off()
pdf(file = "KAPLAN-MEIER CURVE OF GSE31210.pdf",width = 8,height = 6)
compSurv(moic.res         = yau.ntp.pred,
         surv.info        = GSE31210$clin.info,
         convt.time       = "y", # switch to month
         surv.median.line = "hv", # switch to both
         fig.name         = "KAPLAN-MEIER CURVE OF GSE31210") 
dev.off()

markergene=marker.up$templates
write.table(markergene,file = "markergene.txt",sep = "\t",quote = F,row.names = F)









