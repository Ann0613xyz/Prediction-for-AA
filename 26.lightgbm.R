library(limma)
library(caret)
library(pROC)
setwd=("C:\\Users\\anxin\\Desktop\\integrated_model")
set.seed(22)
dat=read.table("C:\\Users\\anxin\\Desktop\\integration\\normalize.txt", header=TRUE, sep="\t", check.names=F)
dat=na.omit(dat)
dat=as.data.frame(t(dat))
diagnosis <- sapply(strsplit(rownames(dat), "_"), `[`, 2)
train=cbind(diagnosis, dat)
gene_symbol <- train[1, -1] 
colnames(train)[-1] <- gene_symbol 
train <- train[-1, ]


outData=read.table("C:\\Users\\anxin\\Desktop\\integration\\testnormalize.txt", header=TRUE, sep="\t", check.names=F)

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
dat.train<-na.omit(dat.train)
dat.test<-na.omit(dat.test)

library(lightgbm)
library(caret)
library(Metrics)
y_train<-as.matrix(dat.train[,1])
x_train <- as.matrix(dat.train[, 2:ncol(dat.train)])

y_train_num <- as.numeric(as.character(dat.train[,1]))
x_train <- as.matrix(dat.train[, 2:ncol(dat.train)])
dtrain <- lgb.Dataset(data = x_train, label = y_train_num)


lgb_grid <- expand.grid(
  num_leaves = c(15, 31, 63),
  learning_rate = c(0.01, 0.05, 0.1),
  stringsAsFactors = FALSE
)


num_rounds <- 100
best_logloss <- Inf
best_params <- NULL


y_train_num <- as.numeric(as.character(dat.train[,1]))
x_train <- as.matrix(dat.train[, 2:ncol(dat.train)])


lgb_grid <- expand.grid(
  num_leaves = c(15, 31, 63),
  learning_rate = c(0.01, 0.05, 0.1),
  stringsAsFactors = FALSE
)


num_rounds <- 100
best_logloss <- Inf
best_params <- NULL


base_params <- list(
  objective = "binary",
  metric = "binary", 
  boosting_type = "gbdt", 
  feature_pre_filter = FALSE, 
  min_data_in_leaf = 1,        
  min_data_in_feature = 1,     
  feature_fraction = 0.9,
  bagging_fraction = 0.8,
  bagging_freq = 5,
  verbose = -1,
  nthread = 2
)

dtrain <- lgb.Dataset(
  data = x_train, 
  label = y_train_num, 
  params = base_params 
)


for (i in 1:nrow(lgb_grid)) {
  current_grid_params <- as.list(lgb_grid[i, ])
  params_cv <- c(current_grid_params, base_params)
  cv_results <- lgb.cv(
    params = params_cv, 
    data = dtrain,
    nrounds = num_rounds,
    nfold = 5,             
    eval_freq = 10
  )
  logloss_mean_vector <- unlist(cv_results$record_evals$valid$binary$mean)
  
  if (length(logloss_mean_vector) > 0) {
    mean_logloss <- min(logloss_mean_vector)
  } else {
    mean_logloss <- Inf
    cat("Warning: This combination did not record a valid LogLoss result!\n")
  }
  
  cat(sprintf("combination %d: num_leaves=%d, lr=%.2f, LogLoss=%.4f\n", 
              i, params_cv$num_leaves, params_cv$learning_rate, mean_logloss))
  
  if (mean_logloss < best_logloss) {
    best_logloss <- mean_logloss
    best_params <- params_cv
  }
}

cat(paste("best LogLoss:", round(best_logloss, 4), "\n"))
print(best_params)

if (!is.null(best_params)) {
  model <- lgb.train(
    params = best_params,
    data = dtrain,
    nrounds = num_rounds
  )
  cat("\nparams got\n")
} else {
  cat("\nfailed to find params\n")}

model <- lgb.train(
  params = best_params,
  data = dtrain,
  nrounds = num_rounds
)

y_test <- as.matrix(dat.test$diagnosis)
x_test <- as.matrix(dat.test[,2:ncol(dat.test)])


y_pred_probs <- predict(model, x_test)
y_pred <- ifelse(y_pred_probs > 0.5, 1, 0)


yt_pred_probs <- predict(model, x_train)
yt_pred <- ifelse(yt_pred_probs > 0.5, 1, 0)


y_train <- as.factor(y_train)
y_test<-as.factor(y_test)

y_pred <- factor(y_pred, levels = levels(y_test))
conf_matrix<-confusionMatrix(y_pred,dat.test$diagnosis)
conf_matrix
conf_matrix_table <- conf_matrix$table


TP <- conf_matrix_table[2, 2]
FP <- conf_matrix_table[1, 2] 
FN <- conf_matrix_table[2, 1] 
TN <- conf_matrix_table[1, 1]
N_total <- sum(conf_matrix_table)

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
  pred = y_pred_probs, 
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


yt_pred <- factor(yt_pred, levels = levels(y_train))
conf_matrix<-confusionMatrix(yt_pred,dat.train$diagnosis)
conf_matrix
conf_matrix_table <- conf_matrix$table

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



train_data_boot <- data.frame(
  pred = yt_pred_probs,
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
          file = "lightGBM_Model_Performance_Metrics_with_CI.csv", 
          row.names = FALSE)

trainData<-dat.train
testData<-dat.test

# ROC_train(text labels should be manually annotated)
predict <- predict(model,as.matrix(trainData[,2:ncol(trainData)]))
pred<-predict
roc.list1 <- roc(trainData$diagnosis, pred)
roc.list1
ci.auc(roc.list1)
g.list1<-ggroc(roc.list1,alpha=1,size=1,legacy.axes=TRUE,color="#1B90CF")
g.list1+theme_update()+
  annotate(geom="segment",x=0,y=0,xend=1,yend=1,linetype=2)+
  annotate("text",x=0.7,y=0.20,
           label="AUC=0.994",colour="#5d6174",size=8)+
  ggtitle("ROC-train")

# ROC_test(text labels should be manually annotated)
predict <- predict(model,as.matrix(testData[,2:ncol(testData)]))
pred<-predict
roc.list2 <- roc(testData$diagnosis, pred)
roc.list2
ci.auc(roc.list2)
g.list1<-ggroc(roc.list2,alpha=1,size=1,legacy.axes=TRUE,color="#1B90CF")
g.list1+theme_update()+
  annotate(geom="segment",x=0,y=0,xend=1,yend=1,linetype=2)+
  annotate("text",x=0.7,y=0.20,
           label="AUC=0.9097",colour="#5d6174",size=8)+
  ggtitle("ROC-test")

#DCA
t1<-dat.train
t2<-dat.test

t1$prob<-yt_pred_probs
t2$prob<-y_pred_probs
library(dcurves)
dca(diagnosis~prob,data=t1,
    as_probability="prob")%>%
  plot(smooth=TRUE)+
  ggplot2::labs(x="Threshold Probability")

dca(diagnosis~prob,data=t2,
    as_probability="prob")%>%
  plot(smooth=TRUE)+
  ggplot2::labs(x="Threshold Probability")
