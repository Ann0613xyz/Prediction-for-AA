library(corrplot)
library(glmnet)
library(caret)
library(nortest)
library(tidyverse)
library(ggpubr)
library(rms)
library(pROC)
library(limma)





set.seed(22)
dat=read.table("C:\\Users\\anxin\\Desktop\\integrated_model\\normalize.txt", header=TRUE, sep="\t", check.names=F)
dat=na.omit(dat)
dat=as.data.frame(t(dat))
diagnosis <- sapply(strsplit(rownames(dat), "_"), `[`, 2)
train=cbind(diagnosis, dat)
gene_symbol <- train[1, -1] 
colnames(train)[-1] <- gene_symbol 
train <- train[-1, ]



outData<-read.table("C:\\Users\\anxin\\Desktop\\integrated_model\\testnormalize.txt", header=TRUE, sep="\t", check.names=F)

rownames(outData) <- outData[, 1]

test <- outData[, -1]
test=as.data.frame(normalizeBetweenArrays(test))
test=as.data.frame(t(test))

test <- cbind(rownames(test), test)
test[, 1] <- sub(".*_", "", test[, 1])

colnames(test)[1] <- "diagnosis"
test$diagnosis=as.factor(test$diagnosis)



train <- na.omit(train)
numeric_vars <- c("KRT83","PPP1R1C","PIRT")
train[numeric_vars] <- lapply(train[numeric_vars], function(x) {
  as.numeric(as.character(x)) 
})
test[numeric_vars] <- lapply(test[numeric_vars], function(x) {
  as.numeric(as.character(x)) 
})
sapply(train[numeric_vars], class)
sapply(test[numeric_vars], class)
for (var in numeric_vars) {
  train[[var]] <- as.numeric(gsub("[^0-9.]", "", train[[var]]))
  test[[var]] <- as.numeric(gsub("[^0-9.]", "", test[[var]]))
}
train <- na.omit(train)
cor_matrix <- cor(train[numeric_vars], use = "complete.obs")
print(cor_matrix)

testRes<-cor.mtest(train[numeric_vars],conf.level=0.95)
testRes
par(mfrow=c(2,3))


corrplot(cor_matrix,method='circle') 
corrplot(
  cor_matrix,
  method='color',
  type='upper',
  add=TRUE,
  tl.pos='n',
  cl.pos='n',
  diag=F,
  p.mat=testRes$p,
  sig.level=c(0.001,0.01,0.05),
  pch.cex=1.5,
  insig='label_sig'
)  

train_numeric <- as.data.frame(train[, numeric_vars])
var_names <- names(train_numeric)
vif_results <- data.frame(Variable = character(), VIF = numeric(), stringsAsFactors = FALSE)

library(car)

for (i in 1:length(var_names)) {
  response <- var_names[i]
  predictors <- var_names[-i]
  formula <- as.formula(paste(response, "~", paste(predictors, collapse = "+")))
  model <- lm(formula, data = train_numeric)
  vif_values <- vif(model)
  vif_results <- rbind(vif_results, data.frame(Variable = response, VIF = vif_values))
}


print(vif_results)
write.csv(vif_results,"C:\\Users\\anxin\\Desktop\\integrated_model\\Feature genes VIF.csv")

train$diagnosis <- as.character(train$diagnosis)
train$diagnosis[train$diagnosis == "con"] <- "0"
train$diagnosis[train$diagnosis == "AA"] <- "1"

train$diagnosis <- factor(train$diagnosis, levels = c("0", "1")) 


test$diagnosis <- as.character(test$diagnosis)
test$diagnosis[test$diagnosis == "con"] <- "0"
test$diagnosis[test$diagnosis == "AA"] <- "1"

test$diagnosis <- factor(test$diagnosis, levels = c("0", "1"))

mydata<-train
attach(mydata)
dd<-datadist(mydata)
options(datadist=dd)

fit0<-lrm(diagnosis~KRT83+PPP1R1C+PIRT,
          data=mydata,x=TRUE,y=TRUE)

fit0 
mydata$diagnosis <- as.factor(mydata$diagnosis)
nom0<-nomogram(fit0,fun=plogis,fun.at=c(.001,.003,.01,.04,.1,.25,.5,.75,.9,.97,.99,.999),
               lp=TRUE,funlabel="diagnosis rate")
plot(nom0)

gd<-predict(fit0,newdata=train,
            se.fit=FALSE,dispersion=NULL,terms=NULL,
            na.action=na.pass)
gd2<-predict(fit0,newdata=test,
             se.fit=FALSE,dispersion=NULL,terms=NULL,
             na.action=na.pass)
library(caret)
library(pROC)
train$predicted <- ifelse(gd > 0.5, 1, 0)
test$predicted <- ifelse(gd2 > 0.5, 1, 0)
confusion_matrix_train <- confusionMatrix(
  factor(train$predicted, levels = c(0, 1), labels = c("0", "1")),
  train$diagnosis
)
print(confusion_matrix_train)
conf_matrix_table <- confusion_matrix_train$table


TP <- conf_matrix_table[2, 2] 
FP <- conf_matrix_table[1, 2] 
FN <- conf_matrix_table[2, 1] 
TN <- conf_matrix_table[1, 1] 
N_total <- sum(conf_matrix_table)
library(Hmisc)

sens_ci <- binconf(TP, TP + FN, alpha = 0.05, method = "wilson")

spec_ci <- binconf(TN, TN + FP, alpha = 0.05, method = "wilson")

prec_ci <- binconf(TP, TP + FP, alpha = 0.05, method = "wilson")

acc_ci <- binconf(TP + TN, N_total, alpha = 0.05, method = "wilson")



accuracy <- acc_ci[1]
sensitivity <- sens_ci[1]
specificity <- spec_ci[1]
precision <- prec_ci[1]
recall <- sens_ci[1]
f1_score <- 2 * precision * recall / (precision + recall)
g_mean <- sqrt(precision * recall)
calc_f1_gmean <- function(data, indices) {
  d <- data[indices, ]
  pred_labels <- ifelse(d[, 1] > 0.5, 1, 0)
  pred_factor <- factor(pred_labels, levels = c(0, 1))
  true_factor <- factor(d[, 2], levels = c(0, 1)) 
  library(caret) 
  cm <- confusionMatrix(pred_factor, true_factor)
  
  cm_table <- cm$table
  
  TP <- ifelse("1" %in% rownames(cm_table) & "1" %in% colnames(cm_table), cm_table["1", "1"], 0)
  FP <- ifelse("1" %in% rownames(cm_table) & "0" %in% colnames(cm_table), cm_table["1", "0"], 0)
  FN <- ifelse("0" %in% rownames(cm_table) & "1" %in% colnames(cm_table), cm_table["0", "1"], 0)
  precision <- ifelse((TP + FP) > 0, TP / (TP + FP), 0)
  recall <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
  f1_score_val <- ifelse((precision + recall) > 0, 2 * precision * recall / (precision + recall), 0)
  g_mean_val <- sqrt(precision * recall)
  
  return(c(f1_score_val, g_mean_val))
}


train_data_boot <- data.frame(
  pred = gd,
  true = train$diagnosis
)


set.seed(22) 
boot_results_train <- boot(data = train_data_boot, statistic = calc_f1_gmean, R = 1000)


f1_ci_train <- boot.ci(boot_results_train, type = "perc", index = 1)
gmean_ci_train <- boot.ci(boot_results_train, type = "perc", index = 2)
train_results <- data.frame(
  Metric = c("Accuracy", "Sensitivity (Recall)", "Specificity", "Precision", "F1 Score", "G-mean"),
  Value = c(accuracy, sensitivity, specificity, precision, f1_score, g_mean),
  Lower_CI = c(acc_ci[2], sens_ci[2], spec_ci[2], prec_ci[2], f1_ci_train$percent[4], gmean_ci_train$percent[4]),
  Upper_CI = c(acc_ci[3], sens_ci[3], spec_ci[3], prec_ci[3], f1_ci_train$percent[5], gmean_ci_train$percent[5])
)

confusion_matrix_test <- confusionMatrix(factor(test$predicted), test$diagnosis)
print(confusion_matrix_test)
conf_matrix_table <- confusion_matrix_test$table


TP <- conf_matrix_table[2, 2]
FP <- conf_matrix_table[1, 2] 
FN <- conf_matrix_table[2, 1] 
TN <- conf_matrix_table[1, 1] 
N_total <- sum(conf_matrix_table)
sens_ci <- binconf(TP, TP + FN, alpha = 0.05, method = "wilson")
spec_ci <- binconf(TN, TN + FP, alpha = 0.05, method = "wilson")
prec_ci <- binconf(TP, TP + FP, alpha = 0.05, method = "wilson")
acc_ci <- binconf(TP + TN, N_total, alpha = 0.05, method = "wilson")


accuracy <- acc_ci[1]
sensitivity <- sens_ci[1]
specificity <- spec_ci[1]
precision <- prec_ci[1]
recall <- sens_ci[1]
f1_score <- 2 * precision * recall / (precision + recall)
g_mean <- sqrt(precision * recall)


test_data_boot <- data.frame(
  pred = gd2, 
  true = test$diagnosis 
)

library(boot)
set.seed(22)
boot_results_test <- boot(data = test_data_boot, statistic = calc_f1_gmean, R = 1000)


f1_ci_test <- boot.ci(boot_results_test, type = "perc", index = 1)
gmean_ci_test <- boot.ci(boot_results_test, type = "perc", index = 2)
test_results <- data.frame(
  Metric = c("Accuracy", "Sensitivity (Recall)", "Specificity", "Precision", "F1 Score", "G-mean"),
  Value = c(accuracy, sensitivity, specificity, precision, f1_score, g_mean),
  Lower_CI = c(acc_ci[2], sens_ci[2], spec_ci[2], prec_ci[2], f1_ci_test$percent[4], gmean_ci_test$percent[4]),
  Upper_CI = c(acc_ci[3], sens_ci[3], spec_ci[3], prec_ci[3], f1_ci_test$percent[5], gmean_ci_test$percent[5])
)
library(tidyverse)
full_results <- merge(train_results, test_results, by="Metric", suffixes = c(".Train", ".Test"))
full_results <- full_results %>%
  mutate(
    `Training Set Value [95% CI]` = sprintf("%.4f [%.4f - %.4f]", Value.Train, Lower_CI.Train, Upper_CI.Train),
    `Validation Set Value [95% CI]` = sprintf("%.4f [%.4f - %.4f]", Value.Test, Lower_CI.Test, Upper_CI.Test)
  ) %>%
  select(Metric, `Training Set Value [95% CI]`, `Validation Set Value [95% CI]`)

print(full_results, row.names = FALSE)
setwd("C:\\Users\\anxin\\Desktop\\integrated_model")
write.csv(full_results, 
          file = "Logit_Model_Performance_Metrics_with_CI.csv", 
          row.names = FALSE)


library(pROC)
library(ggplot2)
# ROC_train(text labels should be manually annotated)
roc.list1<-roc(train$diagnosis,gd)
roc.list1
ci_result1 <- ci.auc(roc.list1)
ci_result1
g.list1<-ggroc(roc.list1,alpha=1,size=1,legacy.axes=TRUE,color="#1B90CF")
g.list1+theme_update()+
  annotate(geom="segment",x=0,y=0,xend=1,yend=1,linetype=2)+
  annotate("text",x=0.7,y=0.20,
           label="AUC=0.936",colour="#5d6174",size=8)+
  ggtitle("ROC-train")
# ROC_test(text labels should be manually annotated)
roc.list2<-roc(test$diagnosis,gd2)
roc.list2
ci_result2 <- ci.auc(roc.list2)
ci_result2
g.list2<-ggroc(roc.list2,alpha=1,size=1,legacy.axes=TRUE,color="#1B90CF")
g.list2+theme_update()+
  annotate(geom="segment",x=0,y=0,xend=1,yend=1,linetype=2)+
  annotate("text",x=0.7,y=0.20,
           label="AUC=0.803",colour="#5d6174",size=8)+
  ggtitle("ROC-test")

#DCA
t1<-train
t2<-test

t1$prob<-gd
t2$prob<-gd2
library(dcurves)
dca(diagnosis~prob,data=t1,
    as_probability="prob")%>%
  plot(smooth=TRUE)+
  ggplot2::labs(x="Threshold Probability")

dca(diagnosis~prob,data=t2,
    as_probability="prob")%>%
  plot(smooth=TRUE)+
  ggplot2::labs(x="Threshold Probability")
