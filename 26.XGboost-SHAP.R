library(corrplot)
library(glmnet)
library(caret)
library(CBCgrps)
library(nortest)
library(tidyverse)
library(ggpubr)
library(rms)
library(pROC)
library(viridis)
library(limma)
library(Hmisc)
library(boot)


setwd=("C:\\Users\\anxin\\Desktop\\integrated_model")


set.seed(22)
dat=read.table("C:\\Users\\anxin\\Desktop\\integrated_model\\normalize.txt", header=TRUE, sep="\t", check.names=F)
dat=na.omit(dat)
dat=as.data.frame(t(dat))
diagnosis <- sapply(strsplit(rownames(dat), "_"), `[`, 2)
train=cbind(diagnosis, dat)
gene_symbol <- train[1, -1] 
colnames(train)[-1] <- gene_symbol 

train <- train[-1, ]

feature_genes <- c("KRT83","PPP1R1C","PIRT")

outData<-read.table("C:\\Users\\anxin\\Desktop\\integrated_model\\testnormalize.txt", header=TRUE, sep="\t", check.names=F)

rownames(outData) <- outData[, 1]

test <- outData[, -1]
test=as.data.frame(normalizeBetweenArrays(test))
test=as.data.frame(t(test))

test <- cbind(rownames(test), test)
test[, 1] <- sub(".*_", "", test[, 1])

colnames(test)[1] <- "diagnosis"
test$diagnosis=as.factor(test$diagnosis)

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


#XGboost
library(data.table)
library(xgboost)
library(Matrix)
library(caTools)


dat.train <- na.omit(dat.train)
x_train <- dat.train[, -1] 
y_train <- dat.train[, 1]

x_train <- as.matrix(data.frame(lapply(x_train, as.numeric)))
y_train <- as.numeric(as.character(y_train))

dtrain <- xgb.DMatrix(data = x_train, label = y_train)


xgb_grid <- expand.grid(
  max_depth = c(2, 3, 4),      
  eta = c(0.01, 0.05, 0.1),  
  gamma = c(0, 0.1),      
  stringsAsFactors = FALSE
)




nrounds <- 100
best_logloss <- Inf
best_params <- NULL


base_params <- list(
  objective = "binary:logistic",
  eval_metric = "logloss",
  min_child_weight = 1, 
  subsample = 0.8,
  colsample_bytree = 0.8, 
  nthread = 2,
  verbose = 0
)

cat("grid search...\n")

for (i in 1:nrow(xgb_grid)) {
  current_grid_params <- as.list(xgb_grid[i, ])
  params_cv <- c(current_grid_params, base_params)
  cv_results <- xgb.cv(
    params = params_cv,
    data = dtrain,
    nrounds = nrounds,
    nfold = 5,
    showsd = FALSE,
    verbose = 0,
    early_stopping_rounds = 10, 
    maximize = FALSE
  )
  
  logloss_mean_vector <- cv_results$evaluation_log$test_logloss_mean

  if (length(logloss_mean_vector) > 0) {
    mean_logloss <- min(logloss_mean_vector)
  } else {
    mean_logloss <- Inf
    cat("Warning: This combination did not record a valid LogLoss result!\n")
  }

  cat(sprintf("combination %d: max_depth=%d, eta=%.2f, gamma=%.1f, LogLoss=%.4f\n", 
              i, params_cv$max_depth, params_cv$eta, params_cv$gamma, mean_logloss))


  if (mean_logloss < best_logloss) {
    best_logloss <- mean_logloss
    best_params <- params_cv
  }
}

cat(paste("Best LogLoss:", round(best_logloss, 4), "\n"))
print(best_params)


if (!is.null(best_params)) {
  xgb_model_final <- xgboost(
    params = best_params,
    data = dtrain,
    nrounds = nrounds
  )
  cat("\nModel has been trained\n")
} else {
  cat("\nModel fail to be trained\n")
}
##validation
train_predictions<-predict(xgb_model_final,newdata=dtrain)
train_predictions1<-ifelse(train_predictions>0.5,1,0)


x_test=dat.test[,-1]
x_test <- data.frame(lapply(x_test, as.numeric))

y_test=dat.test[,1]

dat.test <- na.omit(dat.test)
x_test<-na.omit(x_test)




dat.test[, 1] <- as.numeric(dat.test[, 1])
dtest<-xgb.DMatrix(data=as.matrix(x_test),label=dat.test$diagnosis)
test_predictions<-predict(xgb_model_final,newdata=dtest)
test_predictions1<-ifelse(test_predictions>0.5,1,0)


calc_f1_gmean <- function(data, indices) {

  d <- data[indices, ]
  
  pred_factor <- factor(d[, 1], levels = c(0, 1))
  true_factor <- factor(d[, 2], levels = c(0, 1))
  

  cm <- confusionMatrix(pred_factor, true_factor)
  

  cm_table <- cm$table
  TP <- cm_table[2, 2]
  FP <- cm_table[1, 2]
  FN <- cm_table[2, 1]

  precision <- ifelse((TP + FP) > 0, TP / (TP + FP), 0)
  recall <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)

  f1_score_val <- ifelse((precision + recall) > 0, 2 * precision * recall / (precision + recall), 0)
  g_mean_val <- sqrt(precision * recall)
  

  return(c(f1_score_val, g_mean_val))
}

library(caret)
train_predictions1<-as.factor(train_predictions1)
train[, 1][train[, 1] == "con"] <- 0
train[, 1][train[, 1] == "AA"] <- 1

train$diagnosis <- as.factor(train$diagnosis)

levels(train_predictions1) <- levels(train$diagnosis)
conf_matrix<-confusionMatrix(train_predictions1,train$diagnosis)
conf_matrix
conf_matrix_table <- conf_matrix$table


TP <- conf_matrix_table[2, 2]
FP <- conf_matrix_table[1, 2]
FN <- conf_matrix_table[2, 1]
TN <- conf_matrix_table[1, 1]

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
  pred = train_predictions1,
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


#Power analysis
# install.packages("pwr")
library(pwr)

power_data <- data.frame(diagnosis = dat.train$diagnosis, x_train)

power_data$diagnosis <- factor(power_data$diagnosis, levels = c(0, 1))

glm_model <- glm(diagnosis ~ KRT83 + PPP1R1C + PIRT, data = power_data, family = binomial)

# install.packages("pscl")
library(pscl)
R2_metrics <- pscl::pR2(glm_model)
R2_pseudo <- R2_metrics["McFadden"]
f2_effect_size <- R2_pseudo / (1 - R2_pseudo)
N_train <- nrow(power_data)
u_df <- length(feature_genes)
v_df <- N_train - u_df - 1
alpha_level <- 0.05
power_result <- pwr.f2.test(u = u_df, v = v_df, f2 = f2_effect_size, sig.level = alpha_level)
train_power <- data.frame(
  N = N_train,
  f2_Effect_Size = round(as.numeric(f2_effect_size), 4),
  Power = round(power_result$power, 4),
  Set = "Training"
)




test_predictions1<-as.factor(test_predictions1)
test[, 1] <- as.character(test[, 1])
test[, 1][test[, 1] == "con"] <- 0
test[, 1][test[, 1] == "AA"] <- 1

test$diagnosis <- as.factor(test$diagnosis)
levels(test_predictions1) <- levels(test$diagnosis)
conf_matrix<-confusionMatrix(test_predictions1,test$diagnosis)
conf_matrix
conf_matrix_table <- conf_matrix$table


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
  pred = test_predictions1, 
  true = test$diagnosis 
)

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


power_data_test <- data.frame(diagnosis = dat.test$diagnosis, x_test)

power_data_test$diagnosis <- factor(power_data_test$diagnosis, levels = c(0, 1))

glm_model_test <- glm(diagnosis ~ KRT83 + PPP1R1C + PIRT, data = power_data_test, family = binomial)

R2_metrics_test <- pscl::pR2(glm_model_test)
R2_pseudo_test <- R2_metrics_test["McFadden"]

f2_effect_size_test <- R2_pseudo_test / (1 - R2_pseudo_test)

N_test <- nrow(power_data_test)

u_df_test <- length(feature_genes)

v_df_test <- N_test - u_df_test - 1

alpha_level <- 0.05


power_result_test <- pwr::pwr.f2.test(u = u_df_test, v = v_df_test, f2 = f2_effect_size_test, sig.level = alpha_level)
test_power <- data.frame(
  N = N_test,
  f2_Effect_Size = round(as.numeric(f2_effect_size_test), 4),
  Power = round(power_result_test$power, 4),
  Set = "Validation"
)



full_results <- merge(train_results, test_results, by="Metric", suffixes = c(".Train", ".Test"))
full_results <- full_results %>%
  mutate(
    `Training Set Value [95% CI]` = sprintf("%.4f [%.4f - %.4f]", Value.Train, Lower_CI.Train, Upper_CI.Train),
    `Validation Set Value [95% CI]` = sprintf("%.4f [%.4f - %.4f]", Value.Test, Lower_CI.Test, Upper_CI.Test)
  ) %>%
  select(Metric, `Training Set Value [95% CI]`, `Validation Set Value [95% CI]`)

setwd=("C:\\Users\\anxin\\Deskstop\\integrated_model")
print(full_results, row.names = FALSE)
write.csv(full_results, 
          file = "Model_XGBoost_Performance_Metrics_with_CI.csv", 
          row.names = FALSE)
power_summary <- rbind(train_power, test_power)
print(power_summary, row.names = FALSE)
write.csv(power_summary, 
          file = "Sample_Power_Analysis_Summary.csv", 
          row.names = FALSE)

#ROCtrain
library(pROC)
## ROC_train(text labels should be manually annotated)
train_predictions<-predict(xgb_model_final,newdata=dtrain)
roc.list1<-roc(train$diagnosis, train_predictions)
roc.list1
ci.auc(roc.list1)
g.list1<-ggroc(roc.list1,alpha=1,size=1,legacy.axes=TRUE,color="#1B90CF")
g.list1+theme_update()+
  annotate(geom="segment",x=0,y=0,xend=1,yend=1,linetype=2)+
  annotate("text",x=0.7,y=0.20,
           label="AUC=0.9948(0.9886-1)",colour="#5d6174",size=4)+
  ggtitle("ROC-train")
# ROC_test(text labels should be manually annotated)
test_predictions<-predict(xgb_model_final,newdata=dtest)
roc.list1<-roc(test$diagnosis, test_predictions)
roc.list1
ci.auc(roc.list1)
g.list1<-ggroc(roc.list1,alpha=1,size=1,legacy.axes=TRUE,color="#1B90CF")
g.list1+theme_update()+
  annotate(geom="segment",x=0,y=0,xend=1,yend=1,linetype=2)+
  annotate("text",x=0.7,y=0.20,
           label="AUC=0.9241(0.8423-1)",colour="#5d6174",size=4)+
  ggtitle("ROC-test")
#DCA
t1<-train
t2<-test

t1$prob<-train_predictions
t2$prob<-test_predictions
library(dcurves)
dca(diagnosis~prob,data=t1,
    as_probability="prob")%>%
  plot(smooth=TRUE)+
  ggplot2::labs(x="Threshold Probability")

dca(diagnosis~prob,data=t2,
    as_probability="prob")%>%
  plot(smooth=TRUE)+
  ggplot2::labs(x="Threshold Probability")

#SHAP
library(shapviz)
for (i in 2:ncol(dat.train)) {
  dat.train[[i]] <- as.numeric(dat.train[[i]])
}
shap<-shapviz(xgb_model_final,X_pred=data.matrix(dat.train[,-1]))
control_row_id <- grep("_con", rownames(dat.train))
disease_row_id <- grep("_AA", rownames(dat.train))

#########If you need to randomly selected two samples, run the following commented code
# control_sample_id <- sample(control_row_id, 1)
# disease_sample_id <- sample(disease_row_id, 1)
# cat("Selected Con samples：", rownames(dat.train)[control_sample_id], "\n")
# cat("Selected AA samples：", rownames(dat.train)[disease_sample_id], "\n")
# sv_waterfall(shap,
#              row_id=control_sample_id,
#              fill_colors=c("#E31F1CD7","#246EE3"))
# sv_force(shap,
#          row_id=control_sample_id,
#          max_display=10,
#          fill_colors=c("#E31F1CD7","#246EE3"))
# sv_waterfall(shap,
#              row_id=disease_sample_id,
#              fill_colors=c("#E31F1CD7","#246EE3"))
# sv_force(shap,
#          row_id=disease_sample_id,
#          max_display=10,
#          fill_colors=c("#E31F1CD7","#246EE3"))


# To reproduce the results in the article, run the following uncommented code
control_sample_id <- which(rownames(dat.train) == "GSM1682045_con")
sv_waterfall(shap,
               row_id = control_sample_id,
               fill_colors = c("#F1B129", "#832165"))

sv_force(shap,
           row_id = control_sample_id,
           max_display = 10,
           fill_colors = c("#F1B129", "#832165"))

disease_sample_id <- which(rownames(dat.train) == "GSM2124826_AA")
sv_waterfall(shap,
             row_id = disease_sample_id,
             fill_colors = c("#F1B129", "#832165"))

sv_force(shap,
         row_id = disease_sample_id,
         max_display = 10,
         fill_colors = c("#F1B129", "#832165"))



sv_importance(shap,kind="beeswarm")

sv_dependence(shap,
              v=c("PPP1R1C","PIRT","KRT83"))

shp_i<-shapviz(xgb_model_final,
               X_pred=data.matrix(dat.train[,-1]),
               interactions=TRUE)
sv_interaction(shp_i)+
  theme(axis.text.x=element_text(angle=45,
                                 vjust=1,
                                 hjust=1))

interaction_array <- shp_i$S_inter
feature_names <- colnames(shp_i$X) 
n_features <- length(feature_names)
interaction_results <- list()
k <- 1
for (i in 1:n_features) {
  for (j in 1:n_features) {
    if (i != j) {
      pair_name <- paste0(feature_names[i], "_", feature_names[j])
      shap_values_for_pair <- interaction_array[, i, j]
      mean_abs_interaction_val <- mean(abs(shap_values_for_pair))
      interaction_results[[k]] <- data.frame(
        Feature1 = feature_names[i],
        Feature2 = feature_names[j],
        Interaction_Pair = pair_name,
        Mean_Abs_SHAP_Interaction = mean_abs_interaction_val
      )
      k <- k + 1
    }
  }
}

interaction_df <- do.call(rbind, interaction_results)

interaction_df$Sorted_Pair <- apply(interaction_df[, c("Feature1", "Feature2")], 1, function(x) paste(sort(x), collapse = "_"))

final_interaction_df <- interaction_df %>%
  group_by(Sorted_Pair) %>%
  summarise(
    Mean_Abs_SHAP_Interaction = mean(Mean_Abs_SHAP_Interaction),
    .groups = 'drop'
  ) %>%
  rename(Interaction_Pair = Sorted_Pair) %>%
  arrange(desc(Mean_Abs_SHAP_Interaction))


print(final_interaction_df, row.names = FALSE)

write.csv(final_interaction_df, 
          file = "SHAP_Interaction_Quantification.csv", 
          row.names = FALSE)


saveRDS(xgb_model_final, "xgb_model_final.rds")
feature_genes <- c("KRT83", "PPP1R1C", "PIRT")
saveRDS(feature_genes, "feature_genes.rds")

#Five fold cross validation
set.seed(123)
folds <- sample(1:5, size = nrow(dat.train), replace = TRUE) 

roc_list <- list()

for (i in 1:5) {  
  train_index <- which(folds != i)
  val_index <- which(folds == i)
  
  train_X <- as.matrix(x_train[train_index, ])
  train_y <- dat.train[train_index, "diagnosis"]
  val_X <- as.matrix(x_train[val_index, ])
  val_y <- dat.train[val_index, "diagnosis"]
  
  dtrain_fold <- xgb.DMatrix(data = train_X, label = train_y)
  dval_fold <- xgb.DMatrix(data = val_X, label = val_y)

  xgb_model_fold <- xgboost(params = best_params, data = dtrain_fold, nrounds = nrounds)

  val_predictions <- predict(xgb_model_fold, newdata = dval_fold)

  roc_obj <- roc(val_y, val_predictions)
  roc_list[[i]] <- roc_obj
}


g_list <- list()
for (i in 1:5) { 
  roc_obj <- roc_list[[i]]
  auc_value <- auc(roc_obj)
  

  g <- ggroc(roc_obj, alpha = 1, size = 1, legacy.axes = TRUE, color = "#1B90CF") +
    theme_update() +
    annotate(geom = "segment", x = 0, y = 0, xend = 1, yend = 1, linetype = 2) +
    annotate("text", x = 0.7, y = 0.20, label = paste("AUC=", round(auc_value, 4)), 
             colour = "#5d6174", size = 8) +
    ggtitle(paste("ROC - Fold", i))
  
  g_list[[i]] <- g
}

for (i in 1:5) { 
  print(g_list[[i]])
}
