rm(list = ls())
library(survival)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(CoxBoost)
library(survivalsvm)
library(dplyr)
library(tibble)
library(BART)
library(miscTools)
library(compareC)
library(ggplot2)
library(ggsci)
library(tidyr)
library(ggbreak)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
load("筛选基因表达临床.RDATA")
OSdata=dplyr::filter(OSdata,OS.time > 0)
aaaaaa=sort(unique(OSdata$cohort))
str(OSdata)
for (i in 1:length(aaaaaa)) {
  aaa=aaaaaa[i]
  # 加载数据集
  linshidata <- OSdata %>%
    dplyr::filter(cohort %in% aaa) %>%
    dplyr::select("ID","OS.time","OS",badgene,goodgene) %>%
    na.omit()
  assign(aaa,linshidata)
}


mm <- list(TCGA=TCGA, GSE13213=GSE13213,GSE31210=GSE31210)


# 数据标准化
mm <- lapply(mm,function(x){
  x[,-c(1:3)] <- scale(x[,-c(1:3)])
  return(x)})

result <- data.frame()
# TCGA作为训练集
est_data <- mm$TCGA
# GEO作为验证集
val_data_list <- mm
pre_var <- colnames(est_data)[-c(1:3)]
est_dd <- est_data[, c('OS.time', 'OS', pre_var)]
val_dd_list <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', pre_var)]})

# 设置种子数和节点数，其中节点数可以调整
rf_nodesize <- 5
seed <- 123

########################################################################################################
#1.RSF
## 1-1.RSF
set.seed(seed)
fit <- rfsrc(Surv(OS.time,OS)~., data = est_dd, 
             ntree = 1000, nodesize = rf_nodesize,  #该值建议多调整
             splitrule = 'logrank', 
             importance = T, 
             proximity = T, 
             forest = T, 
             seed = seed)
rs <- lapply(val_dd_list, function(x){cbind(x[, 1:2], RS  = predict(fit, newdata = x)$predicted)})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- 'RSF'
result <- rbind(result, cc)
mean(result$Cindex)
## 1-2.RSF + CoxBoost
set.seed(seed)
fit <- rfsrc(Surv(OS.time, OS)~., data = est_dd, 
             ntree = 1000, nodesize = rf_nodesize,  #该值建议多调整
             splitrule = 'logrank', 
             importance = T, 
             proximity = T, 
             forest = T, 
             seed = seed)
rid <- var.select(object = fit)

rid <- rid$topvars
est_dd2 <- est_data[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', rid)]})
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1, 2)]), 
                            trace=TRUE, start.penalty = 500, parallel = T)

cv.res <- cv.CoxBoost(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1, 2)]), 
                      maxstepno = 500, K = 10, type = "verweij",  penalty = pen$penalty)
fit <- CoxBoost(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1, 2)]), 
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, newdata = x[, -c(1, 2)], newtime = x[, 1],  newstatus = x[, 2], type = "lp")))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + CoxBoost')
result <- rbind(result, cc)

## 1-3.RSF + Enet
set.seed(seed)
fit <- rfsrc(Surv(OS.time, OS)~., data = est_dd, 
             ntree = 1000, nodesize = rf_nodesize, #该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit)

rid <- rid$topvars
est_dd2 <- est_data[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', rid)]})
x1 <- as.matrix(est_dd2[, rid])
x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
for (alpha in seq(0.1, 0.9, 0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2, family = "cox", alpha = alpha, nfolds = 10)
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'link', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('RSF + ', 'Enet', '[α=', alpha, ']')
  result <- rbind(result, cc)
}

## 1-4.RSF + GBM
set.seed(seed)
fit <- rfsrc(Surv(OS.time, OS)~., data = est_dd,
             ntree = 1000, nodesize = rf_nodesize,  #该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit)

rid <- rid$topvars
est_dd2 <- est_data[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', rid)]})
set.seed(seed)
fit <- gbm(formula = Surv(OS.time, OS)~., data = est_dd2, distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 6)

# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~., data = est_dd2, distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 8)
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})
cc <- data.frame(Cindex=sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + ', 'GBM')
result <- rbind(result, cc)

## 1-5.RSF + Lasso
set.seed(seed)
fit <- rfsrc(Surv(OS.time, OS)~., data = est_dd,
             ntree = 1000, nodesize = rf_nodesize, #该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit)

rid <- rid$topvars
est_dd2 <- est_data[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', rid)]})
x1 <- as.matrix(est_dd2[, rid])
x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "binomial", alpha = 1,
                type.measure = "class")
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + ', 'Lasso')
result <- rbind(result, cc)

## 1-6.RSF + plsRcox
set.seed(seed)
fit <- rfsrc(Surv(OS.time, OS)~., data = est_dd,
             ntree = 1000, nodesize = rf_nodesize, #该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit)


rid <- rid$topvars
est_dd2 <- est_data[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', rid)]})
set.seed(seed)
cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[, rid], time = est_dd2$OS.time, status = est_dd2$OS), nt = 10, verbose = FALSE)
fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time, event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))

rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = x[, -c(1, 2)])))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + ', 'plsRcox')
result <- rbind(result, cc)

## 1-7.RSF + Ridge
set.seed(seed)
fit <- rfsrc(Surv(OS.time, OS)~., data = est_dd, 
             ntree = 1000, nodesize = rf_nodesize, #该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit)

rid <- rid$topvars
est_dd2 <- est_data[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', rid)]})
x1 <- as.matrix(est_dd2[, rid])
x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold=10, #例文描述：10-fold cross-validation
                family = "binomial", alpha = 0,
                type.measure = "class")
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + ', 'Ridge')
result <- rbind(result, cc)

## 1-8.RSF + StepCox
set.seed(seed)
fit <- rfsrc(Surv(OS.time,OS)~., data = est_dd,
             ntree = 1000, nodesize = rf_nodesize, #该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit)

rid <- rid$topvars
est_dd2 <- est_data[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', rid)]})
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(OS.time, OS)~., est_dd2), direction = direction)
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS=predict(fit, type = 'risk', newdata = x))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('RSF + ', 'StepCox', '[', direction, ']')
  result <- rbind(result, cc)
}

## 1-9.RSF + SuperPC
set.seed(seed)
fit <- rfsrc(Surv(OS.time,OS)~., data = est_dd,
             ntree = 1000, nodesize = rf_nodesize, ##该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit)

rid <- rid$topvars
est_dd2 <- est_data[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', rid)]})
data <- list(x = t(est_dd2[, -c(1, 2)]), y = est_dd2$OS.time,
             censoring.status = est_dd2$OS, 
             featurenames = colnames(est_dd2)[-c(1, 2)])
set.seed(seed)
fit <- superpc.train(data = data, type = 'survival', s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit, data, n.threshold = 20, #default
                     n.fold = 10, 
                     n.components = 3, 
                     min.features = 5, 
                     max.features = nrow(data$x), 
                     compute.fullcv = TRUE, 
                     compute.preval = TRUE)
rs <- lapply(val_dd_list2, function(w){
  test <- list(x = t(w[, -c(1, 2)]), y = w$OS.time, censoring.status=w$OS, featurenames = colnames(w)[-c(1, 2)])
  ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1, ])], n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[, 1:2], RS = rr)
  return(rr2)
})
cc <- data.frame(Cindex=sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + ', 'SuperPC')
result <- rbind(result, cc)

## 1-10.RSF + survival-SVM
set.seed(seed)
fit <- rfsrc(Surv(OS.time, OS)~., data = est_dd,
             ntree = 1000, nodesize = rf_nodesize, #该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit)

rid <- rid$topvars
est_dd2 <- est_data[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', rid)]})
fit = survivalsvm(Surv(OS.time, OS)~., data= est_dd2, gamma.mu = 1)
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS=as.numeric(predict(fit, x)$predicted))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + ', 'survival-SVM')
result <- rbind(result,cc)


#####################################################################################
#2.Enet
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
for (alpha in seq(0.1, 0.9, 0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2, family = "cox", alpha = alpha, nfolds = 10)
  rs <- lapply(val_dd_list,function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit,type = 'link', newx = as.matrix(x[,-c(1,2)]), s = fit$lambda.min)))})
  cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('Enet', '[α=', alpha, ']')
  result <- rbind(result, cc)
}

####################################################################################
#3.StepCox
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(OS.time,OS)~., est_dd), direction = direction)
  rs <- lapply(val_dd_list,function(x){cbind(x[, 1:2], RS = predict(fit, type = 'risk', newdata = x))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']')
  result <- rbind(result, cc)
}

for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(OS.time, OS)~., est_dd), direction = direction)
  rid <- names(coef(fit))#这里不用卡P值，迭代的结果就是可以纳入的基因
  est_dd2 <- est_data[,c('OS.time', 'OS', rid)]
  val_dd_list2 <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', rid)]})
  set.seed(seed)
  pen <- optimCoxBoostPenalty(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1,2)]),
                              trace=TRUE, start.penalty = 500, parallel = T)
  cv.res <- cv.CoxBoost(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1,2)]),
                        maxstepno = 500, K = 10 , type = "verweij", penalty = pen$penalty)
  fit <- CoxBoost(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1, 2)]),
                  stepno = cv.res$optimal.step, penalty = pen$penalty)
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, newdata = x[, -c(1, 2)], newtime=x[, 1], newstatus=x[,2], type="lp")))})
  cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + CoxBoost')
  result <- rbind(result, cc)
  x1 <- as.matrix(est_dd2[, rid])
  x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
  for (alpha in seq(0.1, 0.9, 0.1)) {
    set.seed(seed)
    fit = cv.glmnet(x1, x2, family = "cox",alpha = alpha, nfolds = 10)
    rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'link', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
    cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
      rownames_to_column('ID')
    cc$Model <- paste0('StepCox', '[', direction, ']', ' + Enet', '[α=', alpha, ']')
    result <- rbind(result, cc)
  }
  set.seed(seed)
  fit <- gbm(formula = Surv(OS.time, OS)~., data = est_dd2, distribution = 'coxph',
             n.trees = 10000,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001,
             cv.folds = 10,n.cores = 6)
  # find index for number trees with minimum CV error
  best <- which.min(fit$cv.error)
  set.seed(seed)
  fit <- gbm(formula = Surv(OS.time, OS)~., data = est_dd2, distribution = 'coxph',
             n.trees = best,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001,
             cv.folds = 10,n.cores = 8)
  rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})
  cc <- data.frame(Cindex=sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + GBM')
  result <- rbind(result, cc)
  x1 <- as.matrix(est_dd2[, rid])
  x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
  set.seed(seed)
  fit = cv.glmnet(x1, x2,
                  nfold=10, #例文描述：10-fold cross-validation
                  family = "binomial", alpha = 1,
                  type.measure = "class")
  rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + Lasso')
  result <- rbind(result, cc)
  set.seed(seed)
  cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[,rid], time = est_dd2$OS.time, status = est_dd2$OS), nt = 10, verbose = FALSE)
  fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time,
                 event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = x[, -c(1,2)])))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + plsRcox')
  result <- rbind(result, cc)
  x1 <- as.matrix(est_dd2[, rid])
  x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
  set.seed(seed)
  fit = cv.glmnet(x1, x2,
                  nfold = 10, #例文描述：10-fold cross-validation
                  family = "binomial", alpha = 0,
                  type.measure = "class")
  rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + Ridge')
  result <- rbind(result, cc)
  set.seed(seed)
  fit <- rfsrc(Surv(OS.time,OS)~., data = est_dd2,
               ntree = 1000, nodesize = rf_nodesize, #该值建议多调整
               splitrule = 'logrank',
               importance = T,
               proximity = T,
               forest = T,
               seed = seed)
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = predict(fit, newdata = x)$predicted)})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + RSF')
  result <- rbind(result, cc)
  data <- list(x = t(est_dd2[, -c(1, 2)]), y = est_dd2$OS.time,
               censoring.status = est_dd2$OS,
               featurenames = colnames(est_dd2)[-c(1,2)])
  set.seed(seed)
  fit <- superpc.train(data = data,type = 'survival', s0.perc = 0.5) #default
  cv.fit <- superpc.cv(fit, data, n.threshold = 20, #default
                       n.fold = 10,
                       n.components = 3,
                       min.features = 5,
                       max.features = nrow(data$x),
                       compute.fullcv = TRUE,
                       compute.preval = TRUE)
  rs <- lapply(val_dd_list2, function(w){
    test <- list(x = t(w[, -c(1,2)]), y = w$OS.time, censoring.status = w$OS, featurenames = colnames(w)[-c(1,2)])
    ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])], n.components = 1)
    rr <- as.numeric(ff$v.pred)
    rr2 <- cbind(w[,1:2], RS = rr)
    return(rr2)
  })
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + SuperPC')
  result <- rbind(result, cc)
  fit = survivalsvm(Surv(OS.time,OS)~., data = est_dd2, gamma.mu = 1)
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x)$predicted))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + survival-SVM')
  result <- rbind(result, cc)
}

#####################################################################
#4.CoxBoost
## 4-1.CoxBoost
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]), 
                            trace = TRUE, start.penalty = 500, parallel = T)

cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]), maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rs <- lapply(val_dd_list, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, newdata = x[, -c(1,2)], newtime = x[,1], newstatus = x[,2], type = "lp")))})
cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost')
result <- rbind(result, cc)

## 4-2.CoxBoost + Enet
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                      maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)`!=0), "id"]
est_dd2 <- est_data[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', rid)]})
x1 <- as.matrix(est_dd2[, rid])
x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
for (alpha in seq(0.1, 0.9, 0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2, family = "cox", alpha = alpha, nfolds = 10)
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'link', newx = as.matrix(x[, -c(1,2)]), s = fit$lambda.min)))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('CoxBoost', ' + Enet', '[α=', alpha, ']')
  result <- rbind(result, cc)
}

## 4-3.CoxBoost + GBM
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)


cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                      maxstepno = 500, K= 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)`!=0), "id"]
est_dd2 <- est_data[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[,c('OS.time', 'OS', rid)]})
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~., data = est_dd2, distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~., data = est_dd2, distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 8)
rs <- lapply(val_dd_list2,function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})
cc <- data.frame(Cindex=sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + ', 'GBM')
result <- rbind(result, cc)

## 4-4.CoxBoost + Lasso
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)


cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                      maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                stepno=cv.res$optimal.step, penalty=pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]
est_dd2 <- est_data[,c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', rid)]})
x1 <- as.matrix(est_dd2[, rid])
x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "binomial", alpha = 1,
                type.measure = "class")
rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1,2)]), s = fit$lambda.min)))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + ', 'Lasso')
result <- rbind(result, cc)

## 4-5.CoxBoost + plsRcox
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)

cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                      maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]
est_dd2 <- est_data[,c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[,c('OS.time', 'OS', rid)]})
set.seed(seed)
cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[,rid], time = est_dd2$OS.time, status = est_dd2$OS), nt = 10, verbose = FALSE)
fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time, event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))
rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type="lp", newdata = x[, -c(1,2)])))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + ', 'plsRcox')
result <- rbind(result, cc)

## 4-6.CoxBoost + Ridge
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)

cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                      maxstepno = 500, K=10, type="verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]
est_dd2 <- est_data[,c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[,c('OS.time', 'OS', rid)]})
x1 <- as.matrix(est_dd2[, rid])
x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold=10, #例文描述：10-fold cross-validation
                family = "binomial", alpha = 0,
                type.measure = "class")
rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1,2)]), s = fit$lambda.min)))})
cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + ', 'Ridge')
result <- rbind(result, cc)

## 4-7.CoxBoost + StepCox
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                      maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]
est_dd2 <- est_data[,c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[,c('OS.time', 'OS', rid)]})
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(OS.time,OS)~., est_dd2), direction = direction)
  rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = predict(fit, type = 'risk', newdata = x))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('CoxBoost + ', 'StepCox', '[', direction, ']')
  result <- rbind(result, cc)
}
## 4-8.CoxBoost + SuperPC
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[,-c(1,2)]),
                      maxstepno = 500, K= 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[,-c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]
est_dd2 <- est_data[,c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[,c('OS.time', 'OS', rid)]})
data <- list(x = t(est_dd2[, -c(1,2)]), y = est_dd2$OS.time, censoring.status = est_dd2$OS,
             featurenames = colnames(est_dd2)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data, type = 'survival', s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit, data, n.threshold = 20, #default
                     n.fold = 10,
                     n.components = 3,
                     min.features = 5,
                     max.features = nrow(data$x),
                     compute.fullcv = TRUE,
                     compute.preval =TRUE)
rs <- lapply(val_dd_list2, function(w){
  test <- list(x=t(w[, -c(1,2)]), y = w$OS.time, censoring.status = w$OS, featurenames = colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])], n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2], RS = rr)
  return(rr2)
})
cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + ', 'SuperPC')
result <- rbind(result, cc)

## 4-9.CoxBoost + survival-SVM
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)

cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                      maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[,-c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]
est_dd2 <- est_data[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', rid)]})
fit = survivalsvm(Surv(OS.time, OS)~., data = est_dd2, gamma.mu = 1)
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x)$predicted))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + ', 'survival-SVM')
result <- rbind(result, cc)





########################################################################################################
#5.plsRcox
set.seed(seed)
cv.plsRcox.res = cv.plsRcox(list(x = est_dd[,pre_var], time = est_dd$OS.time, status = est_dd$OS), nt = 10, verbose = FALSE)
fit <- plsRcox(est_dd[,pre_var], time = est_dd$OS.time, event = est_dd$OS, nt = as.numeric(cv.plsRcox.res[5]))

rs <- lapply(val_dd_list, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit,type = "lp", newdata = x[, -c(1, 2)])))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time,OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('plsRcox')
result <- rbind(result, cc)



########################################################################################################
#6.superpc
data <- list(x = t(est_dd[, -c(1,2)]), y = est_dd$OS.time, censoring.status = est_dd$OS, featurenames = colnames(est_dd)[-c(1, 2)])
set.seed(seed) 
fit <- superpc.train(data = data,type = 'survival', s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit, data, n.threshold = 20, #default
                     n.fold = 10,
                     n.components = 3,
                     min.features = 5,
                     max.features = nrow(data$x),
                     compute.fullcv = TRUE,
                     compute.preval = TRUE)
rs <- lapply(val_dd_list, function(w){
  test <- list(x = t(w[,-c(1,2)]), y = w$OS.time, censoring.status = w$OS, featurenames = colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])], n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2], RS = rr)
  return(rr2)
})
cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('SuperPC')
result <- rbind(result, cc)



########################################################################################################
#7.GBM
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~., data = est_dd, distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time, OS)~., data = est_dd, distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 8)
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('GBM')
result <- rbind(result, cc)



########################################################################################################
#8.survivalsvm

fit = survivalsvm(Surv(OS.time,OS)~., data = est_dd, gamma.mu = 1)
rs <- lapply(val_dd_list, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, x)$predicted))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('survival - SVM')
result <- rbind(result, cc)

########################################################################################################
#9.Ridge
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = glmnet(x1, x2, family = "binomial", alpha = 0, lambda = NULL)
cvfit = cv.glmnet(x1, x2,
                  nfold = 10, #例文描述：10-fold cross-validation
                  family = "binomial",
                  type.measure = "class"
)

rs <- lapply(val_dd_list, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1,2)]), s = cvfit$lambda.min)))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time,OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Ridge')
result <- rbind(result, cc)



########################################################################################################
#当lasso回归，使用family = "binomial"，第一项会有截距，需要运行rid <- rid[-1]去掉，而用family = "cox"的时候就不需要
#10.Lasso
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "cox", alpha = 1,
                type.measure = "class")
rs <- lapply(val_dd_list, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1,2)]), s = fit$lambda.min)))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso')
result <- rbind(result, cc)

## 10.1.Lasso + CoxBoost
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "cox", alpha = 1,
                type.measure = "class")
fit$lambda.min

myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
#rid <- rid[-1]
est_dd2 <- est_data[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', rid)]})
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBoost(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1,2)]),
                      maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, newdata = x[,-c(1,2)], newtime = x[,1], newstatus = x[,2], type = "lp")))})
cc <- data.frame(Cindex=sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + CoxBoost')
result <- rbind(result, cc)

## 10.2.Lasso + GBM
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "cox", alpha = 1,
                type.measure = "class")
fit$lambda.min

myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
#rid <- rid[-1]
est_dd2 <- est_data[,c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[,c('OS.time', 'OS', rid)]})
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~., data = est_dd2, distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~., data = est_dd2, distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 8)
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})
cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + ', 'GBM')
result <- rbind(result, cc)

## 10.3.Lasso + plsRcox
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "cox", alpha = 1,
                type.measure = "class")
fit$lambda.min

myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
#rid <- rid[-1]
est_dd2 <- est_data[,c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[,c('OS.time', 'OS', rid)]})
set.seed(seed)
cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[, rid], time = est_dd2$OS.time, status = est_dd2$OS), nt = 10, verbose = FALSE)
fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time, event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = x[,-c(1,2)])))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + ', 'plsRcox')
result <- rbind(result, cc)

## 10.4.Lasso + RSF
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "cox", alpha = 1,
                type.measure = "class")
fit$lambda.min

myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
#rid<-rid[-1]
est_dd2 <- est_data[,c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', rid)]})
set.seed(seed)
fit <- rfsrc(Surv(OS.time,OS)~., data = est_dd2,
             ntree = 1000, nodesize = rf_nodesize, ##该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = predict(fit, newdata = x)$predicted)})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso', ' + RSF')
result <- rbind(result, cc)

## 10.5.Lasso + stepcox
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "cox", alpha = 1,
                type.measure = "class")
fit$lambda.min

myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
#rid <- rid[-1]
est_dd2 <- est_data[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', rid)]})
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(OS.time,OS)~., est_dd2), direction = direction)
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = predict(fit, type = 'risk', newdata = x))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('Lasso + ', 'StepCox', '[', direction, ']')
  result <- rbind(result, cc)
}

## 10.6.Lasso + superPC
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "cox", alpha = 1,
                type.measure = "class")
fit$lambda.min
myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
#rid <- rid[-1]
est_dd2 <- est_data[,c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', rid)]})
data <- list(x = t(est_dd2[,-c(1,2)]), y = est_dd2$OS.time, censoring.status = est_dd2$OS,
             featurenames = colnames(est_dd2)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival', s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20, #default
                     n.fold = 10,
                     n.components = 3,
                     min.features = 5,
                     max.features = nrow(data$x),
                     compute.fullcv = TRUE,
                     compute.preval = TRUE)
rs <- lapply(val_dd_list2, function(w){
  test <- list(x = t(w[,-c(1,2)]), y = w$OS.time, censoring.status = w$OS, featurenames = colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])], n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[, 1:2], RS = rr)
  return(rr2)
})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + ', 'SuperPC')
result <- rbind(result, cc)

## 10.7.Lasso + survival-SVM
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "cox", alpha = 1,
                type.measure = "class")
fit$lambda.min

myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
#rid <- rid[-1]
est_dd2 <- est_data[,c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[,c('OS.time', 'OS', rid)]})
fit = survivalsvm(Surv(OS.time,OS)~., data = est_dd2, gamma.mu = 1)
rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, x)$predicted))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + ', 'survival-SVM')
result <- rbind(result, cc)


########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################


#将得到的结果赋给result2变量进行操作
result2 <- result
table(result2$Model)
result2$Cindex=as.numeric(result2$Cindex)
###将结果的长数据转换为宽数据
dd2 <- pivot_wider(result2, names_from = 'ID', values_from = 'Cindex') %>% as.data.frame()
str(dd2)

#将C指数定义为数值型
dd2[,-1] <- apply(dd2[,-1], 2, as.numeric)
str(dd2)
#求每个模型的C指数在三个数据集的均值
dd2$All <- apply(dd2[,2:4], 1, mean)
#求每个模型的C指数在GEO验证集的均值
# dd2$GEO <- apply(dd2[,3:4], 1, mean)
###查看每个模型的C指数
head(dd2)


#输出C指数结果
write.table(dd2,"output_C_index.txt", col.names = T, row.names = F, sep = "\t", quote = F)

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
# 根据C指数排序
dd2 <- dd2[order(dd2$All, decreasing = T),]
# 仅绘制GEO验证集的C指数热图

dt <- dd2[, 2:4]
rownames(dt) <- dd2$Model
colnames(dt)
col_ha <- HeatmapAnnotation(which = "col", Cohort = colnames(dt),
                            annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
                            annotation_name_side = "left",
                            #col = list(Cohort=c("Merge"="#00A087B2",
                            #                    "TCGA"="#3C5488B2",
                            #                    "GSE76427" = "blue")),
                            annotation_legend_param = list(Cohort=list(title = "Cohort",
                                                                       title_position = "topleft",
                                                                       title_gp = gpar(fontsize = 12, fontface = "bold"),
                                                                       labels_rot = 0,
                                                                       legend_height = unit(1,"cm"),
                                                                       legend_width = unit(5,"mm"),
                                                                       labels_gp = gpar(fontsize = 9,
                                                                                        fontface = "bold"))
                            )
)
# 行注释
row_ha <- rowAnnotation('Mean Cindex' = anno_barplot(round(rowMeans(dt), 3), bar_width = 1, add_numbers = T,
                                                     labels = c("Mean Cindex"), height = unit(1, "mm"),
                                                     gp = gpar(col = "white", fill = "skyblue1"), numbers_gp = gpar(fontsize = 8),
                                                     axis_param = list(at = c(0, 0.5, 1),
                                                                       labels = c("0", "0.5", "1")),
                                                     width = unit(2.5, "cm")),
                        annotation_name_side = "bottom",
                        annotation_name_gp = gpar(fontsize = 9, fontface = "bold", angle = 90))

# 自定义图形，主要是热图右侧的条形图
cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(
    round(dt[i, j], 2), 
    x, y,
    gp = gpar(
      fontsize = 8
    ))
}

# 画出热图
pdf("ComplexHeatmap.pdf", width = 10, height = 18)
heatmap <- Heatmap(dt,name = " ", #图例标题
                   #参照color_mapping_legend函数设置图例
                   heatmap_legend_param = list(title="",title_position = "topleft", labels_rot = 0,
                                               legend_height = unit(8,"cm"),
                                               legend_width = unit(5,"mm"),
                                               labels_gp = gpar(fontsize = 15, fontface = "bold")),
                   border = TRUE,
                   #column_split = c("Merge","TCGA","GSE76427"),
                   column_gap = unit(3, "mm"),
                   show_column_names = F,
                   show_row_names = T,
                   col = colorRamp2(c(0.4,0.55,0.7), c("#4DBBD5B2", "white", "#E64B35B2")), # 选择颜色
                   column_title ="", # 列标题
                   #row_title ="Intersect Gene",
                   column_title_side = "top", # 列标题位置
                   row_title_side = "left",
                   row_title_rot = 90, # 旋转方向
                   column_title_gp = gpar(fontsize = 12, fontface = "bold",col = "black"), # 颜色，字体，大小
                   #row_title_gp = gpar(fontsize = 15, fontface = "bold",col = "black"),
                   cluster_columns =F,
                   cluster_rows = F,
                   column_order = c(colnames(dt)),
                   show_row_dend = F, # 是否显示聚类树
                   cell_fun = cell_fun,
                   top_annotation = col_ha,
                   right_annotation = row_ha
)
print(heatmap)
dev.off()


set.seed(seed)
fit <- rfsrc(Surv(OS.time,OS)~., data = est_dd, 
             ntree = 1000, nodesize = rf_nodesize,  #该值建议多调整
             splitrule = 'logrank', 
             importance = T, 
             proximity = T, 
             forest = T, 
             seed = seed)
rs <- lapply(val_data_list, function(x){cbind(x[, 1:3], RS  = predict(fit, newdata = x)$predicted)})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')

print(fit)
plot(get.tree(fit,3))
plot(fit)
pdf(file = "SPF model.pdf",width = 8,height = 15)
plot(fit)
dev.off()

## 获得每个样本的 riskscore，进一步进行下游分析######################## Merge

sel=aaaaaa[3]
score_t <- data.frame(get(sel),Score=rs[[sel]]$RS)
#score_t$OS.time=score_t$OS.time/365
library(survival)
library(survminer)
cut <- surv_cutpoint(score_t,'OS.time','OS','Score')
cut
plot(cut)

## 生存分析--Merge TCGA GSE
cat <- surv_categorize(cut)
fit1 <- survfit(Surv(OS.time,OS)~Score,cat)
mytheme <- theme_survminer(font.legend = c(14,"plain", "black"),
                           font.x = c(14,"plain", "black"),
                           font.y = c(14,"plain", "black")) ## 自定义主题

pdf(paste0(sel,".机器学习生存曲线.pdf"),width = 8,height = 8,onefile = F)
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


## ROC
library(timeROC)
col <- c("#0073C2FF","firebrick1","orange") ## 自定义颜色
tt <- timeROC(score_t$OS.time,score_t$OS,score_t$Score,
              cause = 1,weighting = 'marginal',
              times = seq(1,5,1),ROC = T,iid = T)

tt$AUC


pdf(paste0(sel,".时间ROC曲线.pdf"),width = 8,height = 8,onefile = F)
plotAUCcurve(tt,conf.int = T,col = "tomato")
dev.off()

pdf(paste0(sel,".ROC曲线.pdf"),width = 8,height = 8,onefile = F)
plot(tt,time=1,title=FALSE,lwd=1.5,col=col[1])
plot(tt,time=3,col=col[2],add=TRUE,title=FALSE,lwd=1.5)
plot(tt,time=5,col=col[3],add=TRUE,title=FALSE,lwd=1.5)
id <- c(paste0("1-year AUC = ",round(tt$AUC[1],3)),
        paste0("3-year AUC = ",round(tt$AUC[3],3)),
        paste0("5-year AUC = ",round(tt$AUC[5],3))
)
legend("bottomright",id,
       fill=col[1:3],
       bty="o",cex=1,
       border = NA)
abline(0,1,lty=2,lwd=0.5)
dev.off()


save(rs,file = "RS.RDATA")


# 
# 
# #模型对比
# #rm(list = ls())
# load("RS.RDATA")
# 
# 
# genename<-rownames(outTab)
# aaa=as.data.frame(t(outTab))
# aaa$ID=rownames(aaa)
# OSdata=inner_join(clin,aaa,by="ID")
# colnames(OSdata)[2:3]=c("OS.time","OS")
# table(OSdata$cohort)
# OSdata$OS.time=OSdata$OS.time/365
# for (i in 1:length(aaaaaa)) {
#   aaa=aaaaaa[i]
#   # 加载数据集
#   linshidata <- OSdata %>%
#     dplyr::filter(cohort %in% aaa) %>%
#     dplyr::select("ID","OS.time","OS",colnames(OSdata)[11:10489]) %>%
#     na.omit()
#   assign(aaa,linshidata)
# }
# 
# 
# merge= dplyr::select(OSdata,"ID","OS.time","OS",colnames(OSdata)[11:10489])
# merge=na.omit(merge)
# val_dd_list1 <- list(merge=merge,CGGA_301=CGGA_301, CGGA_325=CGGA_325, CGGA_693=CGGA_693,GSE13041=GSE13041,Rembrandt__475=Rembrandt__475,TCGA_702=TCGA_702)
# 
# val_dd_list11 <- lapply(val_data_list, function(x){cbind(x[, c("ID","OS.time","OS")], RS  = predict(fit, newdata = x)$predicted)})
# 
# #PMID36959804
# acgd=c("AHR","MGMT","DACH1","YAP1")
# acgdd=intersect(acgd,colnames(OSdata))
# acgdd
# 
# # 定义要添加的新列的计算公式
# calc_column <- function(x) {
#   x$PMID36959804 = 0.26495*x$AHR +0.30219*x$MGMT-0.08143*x$DACH1+0.07468*x$YAP1
#   return(x)
# }
# 
# val_dd_list1 <- lapply(val_dd_list1,calc_column )
# 
# # PMID36474171
# acgd=c("TUBA1C","RPS4X","KDELR2","SLC40A1")
# acgdd=intersect(acgd,colnames(OSdata))
# acgdd
# 
# calc_column <- function(x) {
#   x$PMID36474171 =  (x$TUBA1C*0.41)+(x$RPS4X*- 0.65)+(x$KDELR2*0.60)
#   return(x)
# }
# 
# val_dd_list1 <- lapply(val_dd_list1,calc_column )
# 
# 
# # PMID36110851
# acgd=c("H19","CYTOR","IGFBP2","KLRC2","C5orf38","CHI3L1")
# acgdd=intersect(acgd,colnames(OSdata))
# acgdd
# 
# calc_column <- function(x) {
#   x$PMID36110851 = 0.2739*x$IGFBP2+0.1406*x$CHI3L1 
#   return(x)
# }
# 
# val_dd_list1 <- lapply(val_dd_list1,calc_column )
# 
# 
# # PMID35982954
# acgd=c("PTX3","BMP2","TNFSF9")
# acgdd=intersect(acgd,colnames(OSdata))
# acgdd
# calc_column <- function(x) {
#   x$PMID35982954 = 0.2576*x$PTX3 +0.2716*x$BMP2 + (-0.2763)*x$TNFSF9 
#   return(x)
# }
# 
# val_dd_list1 <- lapply(val_dd_list1,calc_column )
# 
# 
# # PMID35839032
# acgd=c("RAB42","MSN","TOM1L1","EMP3","APOBEC3C","CRTAC1","FBXO17","HOXA2","SNHG18",
#        "STEAP3","PLAUR","IGFBP2","RARRES1","LOXL1","EN1","FABP5")
# acgdd=intersect(acgd,colnames(OSdata))
# acgdd
# calc_column <- function(x) {
#   x$PMID35839032 = (0.039*x$MSN)+(0.071*x$TOM1L1)+(0.089*x$EMP3)+
#     (0.034*x$APOBEC3C)-(0.035*x$CRTAC1)+(0.013*x$FBXO17) +(0.001*x$HOXA2)+(0.042*x$STEAP3)+(0.077*x$PLAUR)+
#     (0.137*x$IGFBP2)+(0.002*x$RARRES1)+(0.012*x$LOXL1) +(0.099*x$EN1)+(0.017*x$FABP5) 
#   return(x)
# }
# 
# val_dd_list1 <- lapply(val_dd_list1,calc_column )
# 
# 
# # PMID34631528
# acgd=c("COL1A2","COL4A1","ITGB4","MAPK10","PRKCB","PRKCG","RELN","TNC")
# acgdd=intersect(acgd,colnames(OSdata))
# acgdd
# calc_column <- function(x) {
#   x$PMID34631528 = (-0.16036*x$COL1A2)+(0.321935*x$COL4A1)+(0.115143*x$ITGB4)+ 
#     (-0.18388*x$MAPK10)+(-0.20014*x$PRKCB)+(0.191394*x$PRKCG)+ (0.182231*x$RELN)+(0.251663*x$TNC)
#   return(x)
# }
# 
# val_dd_list1 <- lapply(val_dd_list1,calc_column )
# 
# 
# 
# 
# 
# # PMID34589481
# acgd=c("NASP","TIMP1","U2AF2","VEGFA","COL18A1","METTL3","PHLPP2")
# acgdd=intersect(acgd,colnames(OSdata))
# acgdd
# calc_column <- function(x) {
#   x$PMID34589481 = 0.1065*x$NASP + 0.1099*x$TIMP1+ 0.5640*x$U2AF2 + 
#     0.1682*x$VEGFA + (-0.2620)*x$COL18A1 + (-0.1799)*x$METTL3 + (-0.4656)*x$PHLPP2
#   return(x)
# }
# 
# val_dd_list1 <- lapply(val_dd_list1,calc_column )
# 
# 
# # PMID34540896
# acgd=c("HDAC1","RUNX1T1","FKBP3","RBL1","PHF21A")
# acgdd=intersect(acgd,colnames(OSdata))
# acgdd
# calc_column <- function(x) {
#   x$PMID34540896 =0.4013*x$HDAC1-0.3081*x$RUNX1T1-0.4491*x$FKBP3-0.7020*x$PHF21A
#   return(x)
# }
# 
# val_dd_list1 <- lapply(val_dd_list1,calc_column )
# 
# # PMID34123852
# acgd=c("CDK4","HMGB2","WEE1","SMC3","GADD45G")
# acgdd=intersect(acgd,colnames(OSdata))
# acgdd
# calc_column <- function(x) {
#   x$PMID34123852 =0.307*x$CDK4+0.356*x$HMGB2+0.445*x$WEE1+(-0.453*x$SMC3)+(-0.371*x$GADD45G)
#   return(x)
# }
# 
# val_dd_list1 <- lapply(val_dd_list1,calc_column )
# 
# # PMID32296458
# acgd=c("CCL2","IGFBP2","PDPN","SLC12A5")
# acgdd=intersect(acgd,colnames(OSdata))
# acgdd
# calc_column <- function(x) {
#   x$PMID32296458 =(0.2239*x$CCL2) + (0.3375*x$IGFBP2) + (0.1516*x$PDPN) + (0.2276*x$SLC12A5)
#   return(x)
# }
# 
# val_dd_list1 <- lapply(val_dd_list1,calc_column )
# 
# index=c("PMID36959804", "PMID36474171", "PMID36110851", "PMID35982954", "PMID35839032" ,"PMID34631528",
#         "PMID34589481", "PMID34540896", "PMID34123852", "PMID32296458")
# 
# 
# rs_merge=val_dd_list11$merge
# rs_CGGA_301=val_dd_list11$CGGA_301
# rs_CGGA_325=val_dd_list11$CGGA_325
# rs_CGGA_693=val_dd_list11$CGGA_693
# rs_GSE13041=val_dd_list11$GSE13041
# rs_Rembrandt__475=val_dd_list11$Rembrandt__475
# rs_TCGA=val_dd_list11$TCGA_702
# colnames(rs_TCGA)
# 
# rs_merge_cohort=val_dd_list1$merge[,c("ID",index)]
# rs_CGGA_301_cohort=val_dd_list1$CGGA_301[,c("ID",index)]
# rs_CGGA_325_cohort=val_dd_list1$CGGA_325[,c("ID",index)]
# rs_CGGA_693_cohort=val_dd_list1$CGGA_693[,c("ID",index)]
# rs_GSE13041_cohort=val_dd_list1$GSE13041[,c("ID",index)]
# rs_Rembrandt__475_cohort=val_dd_list1$Rembrandt__475[,c("ID",index)]
# rs_TCGA_cohort=val_dd_list1$TCGA_702[,c("ID",index)]
# 
# rs_merge=inner_join(rs_merge,rs_merge_cohort,by = "ID")
# rs_CGGA_301=inner_join(rs_CGGA_301,rs_CGGA_301_cohort,by = "ID")
# rs_CGGA_325=inner_join(rs_CGGA_325,rs_CGGA_325_cohort,by = "ID")
# rs_CGGA_693=inner_join(rs_CGGA_693,rs_CGGA_693_cohort,by = "ID")
# rs_Rembrandt__475=inner_join(rs_Rembrandt__475,rs_Rembrandt__475_cohort,by = "ID")
# rs_TCGA=inner_join(rs_TCGA,rs_TCGA_cohort,by = "ID")
# rs_GSE13041=inner_join(rs_GSE13041,rs_GSE13041_cohort,by = "ID")
# 
# 
# colnames(rs_TCGA)
# list_tcga <- list(RS = rs_TCGA[,c(2:3,4)], PMID36959804  = rs_TCGA[,c(2:3,5)], PMID36474171  = rs_TCGA[,c(2:3,6)],
#                   PMID36110851  = rs_TCGA[,c(2:3,7)], PMID35982954 = rs_TCGA[,c(2:3,8)],
#                   PMID35839032 = rs_TCGA[,c(2:3,9)],  PMID34631528  = rs_TCGA[,c(2:3,10)],
#                   PMID34589481  = rs_TCGA[,c(2:3,11)], PMID34540896  = rs_TCGA[,c(2:3,12)],
#                   PMID34123852 = rs_TCGA[,c(2:3,13)],  PMID32296458 = rs_TCGA[,c(2:3,14)])
#                   
# 
# cc_tcga <- data.frame(Cindex = sapply(list_tcga,function(x){as.numeric(summary(coxph(Surv(OS.time, OS)~., x))$concordance[1])}),
#                       se = sapply(list_tcga,function(x){as.numeric(summary(coxph(Surv(OS.time, OS)~.,x))$concordance[2])})) %>%
#   rownames_to_column('ID')
# 
# 
# colnames(rs_CGGA_301)
# list_CGGA_301 <- list(RS = rs_CGGA_301[,c(2:3,4)], PMID36959804  = rs_CGGA_301[,c(2:3,5)], PMID36474171  = rs_CGGA_301[,c(2:3,6)],
#                   PMID36110851  = rs_CGGA_301[,c(2:3,7)], PMID35982954 = rs_CGGA_301[,c(2:3,8)],
#                   PMID35839032 = rs_CGGA_301[,c(2:3,9)],  PMID34631528  = rs_CGGA_301[,c(2:3,10)],
#                   PMID34589481  = rs_CGGA_301[,c(2:3,11)], PMID34540896  = rs_CGGA_301[,c(2:3,12)],
#                   PMID34123852 = rs_CGGA_301[,c(2:3,13)],  PMID32296458 = rs_CGGA_301[,c(2:3,14)])
# 
# 
# cc_CGGA_301 <- data.frame(Cindex = sapply(list_CGGA_301,function(x){as.numeric(summary(coxph(Surv(OS.time, OS)~., x))$concordance[1])}),
#                       se = sapply(list_CGGA_301,function(x){as.numeric(summary(coxph(Surv(OS.time, OS)~.,x))$concordance[2])})) %>%
#   rownames_to_column('ID')
# 
# colnames(rs_CGGA_325)
# list_CGGA_325 <- list(RS = rs_CGGA_325[,c(2:3,4)], PMID36959804  = rs_CGGA_325[,c(2:3,5)], PMID36474171  = rs_CGGA_325[,c(2:3,6)],
#                   PMID36110851  = rs_CGGA_325[,c(2:3,7)], PMID35982954 = rs_CGGA_325[,c(2:3,8)],
#                   PMID35839032 = rs_CGGA_325[,c(2:3,9)],  PMID34631528  = rs_CGGA_325[,c(2:3,10)],
#                   PMID34589481  = rs_CGGA_325[,c(2:3,11)], PMID34540896  = rs_CGGA_325[,c(2:3,12)],
#                   PMID34123852 = rs_CGGA_325[,c(2:3,13)],  PMID32296458 = rs_CGGA_325[,c(2:3,14)])
# 
# 
# cc_CGGA_325 <- data.frame(Cindex = sapply(list_CGGA_325,function(x){as.numeric(summary(coxph(Surv(OS.time, OS)~., x))$concordance[1])}),
#                       se = sapply(list_CGGA_325,function(x){as.numeric(summary(coxph(Surv(OS.time, OS)~.,x))$concordance[2])})) %>%
#   rownames_to_column('ID')
# 
# colnames(rs_CGGA_693)
# list_CGGA_693 <- list(RS = rs_CGGA_693[,c(2:3,4)], PMID36959804  = rs_CGGA_693[,c(2:3,5)], PMID36474171  = rs_CGGA_693[,c(2:3,6)],
#                   PMID36110851  = rs_CGGA_693[,c(2:3,7)], PMID35982954 = rs_CGGA_693[,c(2:3,8)],
#                   PMID35839032 = rs_CGGA_693[,c(2:3,9)],  PMID34631528  = rs_CGGA_693[,c(2:3,10)],
#                   PMID34589481  = rs_CGGA_693[,c(2:3,11)], PMID34540896  = rs_CGGA_693[,c(2:3,12)],
#                   PMID34123852 = rs_CGGA_693[,c(2:3,13)],  PMID32296458 = rs_CGGA_693[,c(2:3,14)])
# 
# 
# cc_CGGA_693 <- data.frame(Cindex = sapply(list_CGGA_693,function(x){as.numeric(summary(coxph(Surv(OS.time, OS)~., x))$concordance[1])}),
#                       se = sapply(list_CGGA_693,function(x){as.numeric(summary(coxph(Surv(OS.time, OS)~.,x))$concordance[2])})) %>%
#   rownames_to_column('ID')
# 
# colnames(rs_GSE13041)
# list_GSE13041 <- list(RS = rs_GSE13041[,c(2:3,4)], PMID36959804  = rs_GSE13041[,c(2:3,5)], PMID36474171  = rs_GSE13041[,c(2:3,6)],
#                   PMID36110851  = rs_GSE13041[,c(2:3,7)], PMID35982954 = rs_GSE13041[,c(2:3,8)],
#                   PMID35839032 = rs_GSE13041[,c(2:3,9)],  PMID34631528  = rs_GSE13041[,c(2:3,10)],
#                   PMID34589481  = rs_GSE13041[,c(2:3,11)], PMID34540896  = rs_GSE13041[,c(2:3,12)],
#                   PMID34123852 = rs_GSE13041[,c(2:3,13)],  PMID32296458 = rs_GSE13041[,c(2:3,14)])
# 
# 
# cc_GSE13041 <- data.frame(Cindex = sapply(list_GSE13041,function(x){as.numeric(summary(coxph(Surv(OS.time, OS)~., x))$concordance[1])}),
#                       se = sapply(list_GSE13041,function(x){as.numeric(summary(coxph(Surv(OS.time, OS)~.,x))$concordance[2])})) %>%
#   rownames_to_column('ID')
# 
# colnames(rs_Rembrandt__475)
# list_Rembrandt__475 <- list(RS = rs_Rembrandt__475[,c(2:3,4)], PMID36959804  = rs_Rembrandt__475[,c(2:3,5)], PMID36474171  = rs_Rembrandt__475[,c(2:3,6)],
#                   PMID36110851  = rs_Rembrandt__475[,c(2:3,7)], PMID35982954 = rs_Rembrandt__475[,c(2:3,8)],
#                   PMID35839032 = rs_Rembrandt__475[,c(2:3,9)],  PMID34631528  = rs_Rembrandt__475[,c(2:3,10)],
#                   PMID34589481  = rs_Rembrandt__475[,c(2:3,11)], PMID34540896  = rs_Rembrandt__475[,c(2:3,12)],
#                   PMID34123852 = rs_Rembrandt__475[,c(2:3,13)],  PMID32296458 = rs_Rembrandt__475[,c(2:3,14)])
# 
# 
# cc_Rembrandt__475 <- data.frame(Cindex = sapply(list_Rembrandt__475,function(x){as.numeric(summary(coxph(Surv(OS.time, OS)~., x))$concordance[1])}),
#                       se = sapply(list_Rembrandt__475,function(x){as.numeric(summary(coxph(Surv(OS.time, OS)~.,x))$concordance[2])})) %>%
#   rownames_to_column('ID')
# 
# colnames(rs_merge)
# list_merge <- list(RS = rs_merge[,c(2:3,4)], PMID36959804  = rs_merge[,c(2:3,5)], PMID36474171  = rs_merge[,c(2:3,6)],
#                   PMID36110851  = rs_merge[,c(2:3,7)], PMID35982954 = rs_merge[,c(2:3,8)],
#                   PMID35839032 = rs_merge[,c(2:3,9)],  PMID34631528  = rs_merge[,c(2:3,10)],
#                   PMID34589481  = rs_merge[,c(2:3,11)], PMID34540896  = rs_merge[,c(2:3,12)],
#                   PMID34123852 = rs_merge[,c(2:3,13)],  PMID32296458 = rs_merge[,c(2:3,14)])
# 
# 
# cc_merge <- data.frame(Cindex = sapply(list_merge,function(x){as.numeric(summary(coxph(Surv(OS.time, OS)~., x))$concordance[1])}),
#                       se = sapply(list_merge,function(x){as.numeric(summary(coxph(Surv(OS.time, OS)~.,x))$concordance[2])})) %>%
#   rownames_to_column('ID')
# 
# 
# 
# 
# 
# 
# 
# ###第二步，计算其他模型C指数与RS差异的显著性
# reee=c("rs_TCGA","rs_CGGA_301","rs_merge","rs_Rembrandt__475","rs_GSE13041","rs_CGGA_693","rs_CGGA_325")
# reee[2]
# rt <- rs_TCGA
# tcga_compareC_p <- data.frame(Var = colnames(rt[, 4:14]), pval = c(1:length(colnames(rt[, 4:14]))))
# for (i in colnames(rt[, 4:14])) {
#   p <- compareC(rt$OS.time, rt$OS, rt$RS, rt[,i])$pval
#   tcga_compareC_p[which(tcga_compareC_p$Var == i), 2] <- p
# }
# 
# 
# 
# for (a in 1:length(reee)) {
#   rt <- get(reee[a])
#   tcga_compareC_p <- data.frame(Var = colnames(rt[, 4:14]), pval = c(1:length(colnames(rt[, 4:14]))))
#   for (i in colnames(rt[, 4:14])) {
#     p <- compareC(rt$OS.time, rt$OS, rt$RS, rt[,i])$pval
#     tcga_compareC_p[which(tcga_compareC_p$Var == i), 2] <- p
#   }
#   tcga_compareC_p$pval=as.numeric(tcga_compareC_p$pval)
#   tcga_compareC_p[1,2]=1
#   write.table(tcga_compareC_p,paste0("output_",reee[a],"_Cindex_p.txt"), col.names = T, row.names = F, sep = "\t", quote = F)
#   assign(paste0(reee[a],"_compareC_p"),tcga_compareC_p)
# }
# 
# 
# ###第三步，将C指数，置信区间，差异显著性整合
# 
# #合并两个数据
# all_tcga <- data.frame(cc_tcga, p = rs_TCGA_compareC_p[, 2])
# all_CGGA_301 <- data.frame(cc_CGGA_301, p = rs_CGGA_301_compareC_p[, 2])
# all_merge <- data.frame(cc_merge, p = rs_merge_compareC_p[, 2])
# all_CGGA_325 <- data.frame(cc_CGGA_325, p = rs_CGGA_325_compareC_p[, 2])
# all_CGGA_693 <- data.frame(cc_CGGA_693, p = rs_CGGA_693_compareC_p[, 2])
# all_Rembrandt__475 <- data.frame(cc_Rembrandt__475, p = rs_Rembrandt__475_compareC_p[, 2])
# all_GSE13041 <- data.frame(cc_GSE13041, p = rs_GSE13041_compareC_p[, 2])
# 
# 
# 
# 
# ###第四步，画图
# indexx=c("all_tcga","all_CGGA_301","all_merge","all_CGGA_325","all_CGGA_693","all_Rembrandt__475","all_GSE13041")
# for (i in 1:length(indexx)) {
#   dd <- get(indexx[i])
#   #将p值转换为*号
#   dd$ll <- ifelse(dd$p < 0.0001, '****', ifelse(dd$p < 0.001, '***', ifelse(dd$p < 0.01, '**', ifelse(dd$p < 0.05, '*', ''))))
#   rownames(dd) <- NULL
#   library(tidyverse)
#   title=str_replace_all(indexx,"all_","")
#   ggplot(dd, aes(Cindex, reorder(ID, Cindex))) +
#     geom_errorbarh(aes(xmax = Cindex + se, xmin = Cindex - se), color = "black", height = 0, size = 0.7) + # 展示置信区间
#     geom_point(size = 4, shape = 21, fill = pal_nejm()(10)[1]) + # 绘制点图
#     ylab(NULL) + xlab(NULL) +
#     labs(title =title[i]) +
#     geom_vline(xintercept = 0.6, linetype = 'dashed', size = 0.5, color = 'grey50') + # 加纵向虚线
#     theme_bw(base_rect_size = 1) +
#     theme(panel.grid = element_blank(),
#           axis.text.y = element_text(size = 12),
#           axis.text.x = element_text(size = 12),
#           axis.title = element_text(size = 13),
#           plot.title = element_text(hjust = 0.5, size = 15),
#           legend.position = 'none',
#           strip.text = element_text(size = 14)) + 
#     geom_text(aes(x = 0.89, y = ID, label = ll), color = 'black', size = 3, vjust = 0.76) + # 添加*号文本
#     scale_x_continuous(breaks = c(0.5,0.7,0.9), limits = c(0.4, 0.94)) # 设置x轴刻度
#   ggsave(paste0(title[i],"_model_compare.pdf"), width = 4, height = 8)
#   
# }
# 
# save(rs_TCGA,rs_CGGA_301,rs_merge,rs_Rembrandt__475,rs_GSE13041,rs_CGGA_693,rs_CGGA_325,file = "All cohorts results.RDATA")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
