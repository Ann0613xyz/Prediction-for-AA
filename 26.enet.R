library(limma)
library(caret)
library(pROC)
library(glmnet)
setwd("C:\\Users\\anxin\\Desktop\\integrated_model")
set.seed(22)
expr="normalize.txt"
testexpr="testnormalize.txt"
dat=read.table(expr, header=TRUE, sep="\t", check.names=F)
dat=na.omit(dat)
dat=as.data.frame(t(dat))
diagnosis <- sapply(strsplit(rownames(dat), "_"), `[`, 2)
train=cbind(diagnosis, dat)
gene_symbol <- train[1, -1] 
colnames(train)[-1] <- gene_symbol 
train <- train[-1, ]


outData=read.table(testexpr, header=TRUE, sep="\t", check.names=F)

rownames(outData) <- outData[, 1]

test <- outData[, -1]
test=as.data.frame(normalizeBetweenArrays(test))
test=as.data.frame(t(test))

test <- cbind(rownames(test), test)
test[, 1] <- sub(".*_", "", test[, 1])

colnames(test)[1] <- "diagnosis"
train[, 1][train[, 1] == "con"] <- 0
train[, 1][train[, 1] == "AA"] <- 1
test[, 1][test[, 1] == "con"] <- 0
test[, 1][test[, 1] == "AA"] <- 1
train$diagnosis<-as.factor(train$diagnosis)
test$diagnosis<-as.factor(test$diagnosis)

feature_genes <- c("KRT83","PPP1R1C","PIRT")
dat.train=train[,feature_genes]
dat.train <- cbind(rownames(dat.train), dat.train)
colnames(dat.train)[1] <- "diagnosis"
dat.train[, 1] <- sub(".*_", "", dat.train[, 1])
dat.train[, 1][dat.train[, 1] == "con"] <- 0
dat.train[, 1][dat.train[, 1] == "AA"] <- 1

dat.test=test[,feature_genes]
dat.test <- cbind(rownames(dat.test), dat.test)
colnames(dat.test)[1] <- "diagnosis"
dat.test[, 1] <- sub(".*_", "", dat.test[, 1])
dat.test[, 1][dat.test[, 1] == "con"] <- 0
dat.test[, 1][dat.test[, 1] == "AA"] <- 1
dat.train[, 1] <- as.factor(dat.train[, 1])
dat.test[, 1] <- as.factor(dat.test[, 1])

x_train <- as.matrix(dat.train[, 2:ncol(dat.train)])
y_train<-as.matrix(dat.train$diagnosis)
x_test <- as.matrix(dat.test[, 2:ncol(dat.test)])
y_test<-as.matrix(dat.test$diagnosis)

alpha_values <- seq(0, 1, by = 0.1)
best_alpha <- NULL
best_lambda <- NULL
best_auc <- 0

for (alpha in alpha_values) {
  cv_fit <- cv.glmnet(x_train, y_train, family = "binomial", alpha = alpha, nfolds = 10)
  train_pred_probs <- predict(cv_fit, newx = x_train, s = cv_fit$lambda.min, type = "response")
  

  library(pROC)
  auc_value <- auc(y_train, as.numeric(train_pred_probs))
  

  if (auc_value > best_auc) {
    best_auc <- auc_value
    best_alpha <- alpha
    best_lambda <- cv_fit$lambda.min
  }
}


final_model <- glmnet(x_train, y_train, family = "binomial", alpha = best_alpha, lambda = best_lambda)



test_pred_probs <- predict(final_model, newx = x_test, s = best_lambda, type = "response")


test_pred_labels <- ifelse(test_pred_probs > 0.5, 1, 0)


confusion_matrix <- table(test_pred_labels, y_test)
confusion_matrix


TP <- confusion_matrix[2, 2] 
FP <- confusion_matrix[1, 2]
FN <- confusion_matrix[2, 1] 
TN <- confusion_matrix[1, 1] 
N_total <- sum(confusion_matrix)
library(Hmisc)
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
  pred = test_pred_probs,
  true = test$diagnosis   
)

# Bootstrap
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



train_pred_probs <- predict(final_model, newx = x_train, s = best_lambda, type = "response")

train_pred_labels <- ifelse(train_pred_probs > 0.5, 1, 0)

confusion_matrix <- table(train_pred_labels, y_train)
confusion_matrix

TP <- confusion_matrix[2, 2] 
FP <- confusion_matrix[1, 2]
FN <- confusion_matrix[2, 1] 
TN <- confusion_matrix[1, 1] 
N_total <- sum(confusion_matrix)


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



train_data_boot <- data.frame(
  pred = train_pred_probs, 
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
          file = "Enet_Model_Performance_Metrics_with_CI.csv", 
          row.names = FALSE)



# ROC_train(text labels should be manually annotated)
trainData<-dat.train
testData<-dat.test
predict <- predict(final_model,as.matrix(trainData[,2:ncol(trainData)]))
pred<-as.numeric(predict)
roc.list1 <- roc(trainData$diagnosis, pred)
roc.list1
ci.auc(roc.list1)
g.list1<-ggroc(roc.list1,alpha=1,size=1,legacy.axes=TRUE,color="#1B90CF")
g.list1+theme_update()+
  annotate(geom="segment",x=0,y=0,xend=1,yend=1,linetype=2)+
  annotate("text",x=0.7,y=0.20,
           label="AUC=0.936",colour="#5d6174",size=8)+
  ggtitle("ROC-train")

# ROC_test(text labels should be manually annotated)
predict <- predict(final_model,as.matrix(testData[,2:ncol(testData)]))
pred<-as.numeric(predict)
roc.list2 <- roc(testData$diagnosis, pred)
roc.list2
ci.auc(roc.list2)
g.list1<-ggroc(roc.list2,alpha=1,size=1,legacy.axes=TRUE,color="#1B90CF")
g.list1+theme_update()+
  annotate(geom="segment",x=0,y=0,xend=1,yend=1,linetype=2)+
  annotate("text",x=0.7,y=0.20,
           label="AUC=0.803",colour="#5d6174",size=8)+
  ggtitle("ROC-test")

#DCA
t1<-dat.train
t2<-dat.test

t1$prob<-train_pred_probs
t2$prob<-test_pred_probs
library(dcurves)
dca(diagnosis~prob,data=t1,
    as_probability="prob")%>%
  plot(smooth=FALSE)+
  ggplot2::labs(x="Threshold Probability")

dca(diagnosis~prob,data=t2,
    as_probability="prob")%>%
  plot(smooth=FALSE)+
  ggplot2::labs(x="Threshold Probability")
