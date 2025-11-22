library(limma)
library(caret)
library(pROC)
library(glmnet)
dat=read.table("C:\\Users\\anxin\\Desktop\\integrated_model\\normalize.txt", header=TRUE, sep="\t", check.names=F)
dat=na.omit(dat)
dat=as.data.frame(t(dat))
diagnosis <- sapply(strsplit(rownames(dat), "_"), `[`, 2)
train=cbind(diagnosis, dat)
gene_symbol <- train[1, -1] 
colnames(train)[-1] <- gene_symbol 
train <- train[-1, ]


outData=read.table("C:\\Users\\anxin\\Desktop\\integrated_model\\testnormalize.txt", header=TRUE, sep="\t", check.names=F)

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

setwd=("C:\\Users\\anxin\\Desktop\\integrated_model\\model_comparism")

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
  objective = "binary:Logistic",
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
  
  cat(sprintf("Combination %d: max_depth=%d, eta=%.2f, gamma=%.1f, LogLoss=%.4f\n", 
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
dat.test[, 1] <- as.character(dat.test[, 1])
dat.test[, 1] <- as.numeric(dat.test[, 1])
dtest<-xgb.DMatrix(data=as.matrix(x_test),label=dat.test$diagnosis)
test_predictions<-predict(xgb_model_final,newdata=dtest)
test_predictions1<-ifelse(test_predictions>0.5,1,0)
#ROC
library(pROC)
train_predictions<-predict(xgb_model_final,newdata=dtrain)
roc.list_xgbtr<-roc(train$diagnosis, train_predictions)

test_predictions<-predict(xgb_model_final,newdata=dtest)
roc.list_xgbt<-roc(test$diagnosis, test_predictions)

#KNN
y_train<-as.matrix(data.frame(dat.train[,1]))
x_train <- as.matrix(data.frame(dat.train[, 2:ncol(dat.train)]))
y_test<-as.matrix(data.frame(dat.test[,1]))
x_test <- as.matrix(data.frame(dat.test[, 2:ncol(dat.test)]))
set.seed(225)
library(kknn)
knn_train<-knn3(as.matrix(x_train),as.factor(dat.train$diagnosis),k=5)
knn_test<-predict(knn_train,x_test)
knn_test<-apply(knn_test,1,which.max)
knn_train<-predict(knn_train,x_train,type="prob")
knn_train<-apply(knn_train,1,which.max)
roc.list_knntr<-roc(dat.train$diagnosis,as.numeric(knn_train))

roc.list_knnt<-roc(dat.test$diagnosis,as.numeric(knn_test))

#lightGBM
library(lightgbm)
library(caret)
library(Metrics)
y_train<-as.matrix(dat.train[,1])
x_train <- as.matrix(dat.train[, 2:ncol(dat.train)])
dtrain <- lgb.Dataset(data = x_train, label = y_train)
params <- list(
  objective = "binary",
  metric = "binary_logloss",
  learning_rate = 0.1,
  num_leaves = 31,
  feature_fraction = 0.9,
  bagging_fraction = 0.8,
  bagging_freq = 5,
  verbose = -1
)

num_rounds <- 100
model <- lgb.train(
  params = params,
  data = dtrain,
  nrounds = num_rounds
)


y_test <- as.matrix(dat.test$diagnosis)
x_test <- as.matrix(dat.test[,2:ncol(dat.test)])

yt_pred_probs <- predict(model, x_train)
yt_pred <- ifelse(yt_pred_probs > 0.5, 1, 0)
y_pred_probs <- predict(model, x_test)
predict <- predict(model,as.matrix(dat.train[,2:ncol(dat.train)]))
pred<-predict
roc.list_gbmtr <- roc(dat.train$diagnosis, pred)

predict <- predict(model,as.matrix(dat.test[,2:ncol(dat.test)]))
pred<-predict
roc.list_gbmt <- roc(dat.test$diagnosis, pred)

#Enet

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

trainData<-dat.train
testData<-dat.test
predict <- predict(final_model,as.matrix(trainData[,2:ncol(trainData)]))
pred<-as.numeric(predict)
roc.list_enettr <- roc(trainData$diagnosis, pred)

predict <- predict(final_model,as.matrix(testData[,2:ncol(testData)]))
pred<-as.numeric(predict)
roc.list_enett <- roc(testData$diagnosis, pred)

#LR
library(rms)
for (var in feature_genes) {
  train[[var]] <- as.numeric(gsub("[^0-9.]", "", train[[var]]))
  test[[var]] <- as.numeric(gsub("[^0-9.]", "", test[[var]]))
}
mydata<-train
attach(mydata)
dd<-datadist(mydata)
options(datadist=dd)

fit0<-lrm(diagnosis~KRT83+PPP1R1C+PIRT,
          data=mydata,x=TRUE,y=TRUE)
gd<-predict(fit0,newdata=train,
            se.fit=FALSE,dispersion=NULL,terms=NULL,
            na.action=na.pass)
gd2<-predict(fit0,newdata=test,
             se.fit=FALSE,dispersion=NULL,terms=NULL,
             na.action=na.pass)
roc.list_logtr<-roc(train$diagnosis,gd)
roc.list_logt<-roc(test$diagnosis,gd2)

#learning curve
# #comparism

library(data.table)
library(ggplot2)
library(reshape2)
library(caret) 

x_train_temp <- data.frame(lapply(dat.train[, -1], as.numeric))
x_test_temp <- data.frame(lapply(dat.test[, -1], as.numeric))

x_train_full <- as.matrix(x_train_temp)
x_test_full <- as.matrix(x_test_temp)
y_train_full <- as.numeric(as.character(dat.train[, 1]))
y_test_full <- as.numeric(as.character(dat.test[, 1]))

lgb_params <- list(objective = "binary", metric = "binary_logloss", learning_rate = 0.1, num_leaves = 31, verbose = -1)
xgb_params <- list(objective = "binary:Logistic", eval_metric = "logloss", eta = 0.1, max_depth = 3, verbose = 0)
nrounds_fixed <- 100
k_knn <- 5 

enet_alpha <- best_alpha
enet_lambda <- best_lambda

train_sizes_pct <- seq(0.1, 1.0, by = 0.1)
train_sizes_n <- round(train_sizes_pct * nrow(dat.train))
train_sizes_n <- unique(train_sizes_n[train_sizes_n > 0])

learning_curve_data <- data.frame()

for (n_size in train_sizes_n) {

  set.seed(42) 
  sample_indices <- sample(1:nrow(dat.train), size = n_size)
  
  x_sub_train <- x_train_full[sample_indices, ]
  y_sub_train <- y_train_full[sample_indices]

  dtrain_xgb <- xgb.DMatrix(data = x_sub_train, label = y_sub_train)
  dtrain_lgb <- lgb.Dataset(data = x_sub_train, label = y_sub_train)

  
  # 1. XGBoost
  xgb_model <- xgboost(params = xgb_params, data = dtrain_xgb, nrounds = nrounds_fixed)
  pred_xgb <- predict(xgb_model, newdata = xgb.DMatrix(data = x_test_full, label = y_test_full))
  
  # 2. lightGBM
  lgb_model <- lgb.train(params = lgb_params, data = dtrain_lgb, nrounds = nrounds_fixed)
  pred_lgb <- predict(lgb_model, x_test_full)
  
  # 3. KNN
  knn_model <- knn3(x_sub_train, as.factor(y_sub_train), k = k_knn)
  pred_knn <- predict(knn_model, x_test_full, type = "prob")[, 2] 
  
  # 4. Enet
  enet_model <- glmnet(x_sub_train, y_sub_train, family = "binomial", alpha = enet_alpha, lambda = enet_lambda)
  pred_enet <- predict(enet_model, newx = x_test_full, type = "response")
  
  # 5. LR Regression
  df_sub_train <- data.frame(diagnosis = y_sub_train, x_sub_train)
  names(df_sub_train)[-1] <- feature_genes
  
  glm_model <- glm(diagnosis ~ ., data = df_sub_train, family = binomial(link = "logit"))

  df_test_full <- data.frame(x_test_full)
  names(df_test_full) <- feature_genes
  pred_glm <- predict(glm_model, newdata = df_test_full, type = "response")
  
  # LogLoss 
  LogLoss <- function(y_true, y_prob) {
    epsilon <- 1e-15
    y_prob <- pmax(pmin(y_prob, 1 - epsilon), epsilon)
    -mean(y_true * log(y_prob) + (1 - y_true) * log(1 - y_prob))
  }
  
  results <- data.frame(
    N_Train = n_size,
    Model = c("XGBoost", "lightGBM", "KNN", "Enet", "LR"),
    LogLoss = c(
      LogLoss(y_test_full, pred_xgb),
      LogLoss(y_test_full, pred_lgb),
      LogLoss(y_test_full, pred_knn),
      LogLoss(y_test_full, pred_enet),
      LogLoss(y_test_full, pred_glm)
    )
  )
  
  learning_curve_data <- rbind(learning_curve_data, results)
}

colors_lc <- c(XGBoost = "#FF4A46", KNN = "#00BDcd", lightGBM = "#008941", Enet = "#993399", LR = "gold")

g_lc <- ggplot(learning_curve_data, aes(x = N_Train, y = LogLoss, color = Model)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = colors_lc) +
  theme_minimal() +
  labs(
    title = "Learning curves (based on train set)",
    x = "sample size of train set (N)",
    y = "LogLoss for test set",
    color = "model"
  ) +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))

print(g_lc)

#ROC-train
library(pROC)
library(ggplot2)
g.list_xgbtr<-ggroc(roc.list_xgbtr,alpha=1,size=1,legacy.axes=TRUE,color="#FF4A46")
g.list_knntr<-ggroc(roc.list_knntr,alpha=1,size=1,legacy.axes=TRUE,color="#00BDcd")
g.list_gbmtr<-ggroc(roc.list_gbmtr,alpha=1,size=1,legacy.axes=TRUE,color="#008941")
g.list_enettr<-ggroc(roc.list_enettr,alpha=1,size=1,legacy.axes=TRUE,color="#993399")
g.list_logtr<-ggroc(roc.list_logtr,alpha=1,size=1,legacy.axes=TRUE,color="gold")
roc_list <- list(
  Enet = roc.list_enettr,
  lightGBM = roc.list_gbmtr,
  KNN = roc.list_knntr,
  XGBoost = roc.list_xgbtr,
  LR=roc.list_logtr
)

colors <- c(XGBoost = "#FF4A46", KNN = "#00BDcd", lightGBM = "#008941", Enet = "#993399",LR="gold")
roc_data <- do.call(rbind, lapply(names(roc_list), function(name) {
  roc_obj <- roc_list[[name]]
  data.frame(
    FPR = 1 - roc_obj$specificities,
    TPR = roc_obj$sensitivities,
    Model = name
  )
}))


g <- ggplot(roc_data, aes(x = FPR, y = TPR, color = Model)) +
  geom_line(size = 1, alpha = 1) +
  scale_color_manual(values = colors) + 
  theme_minimal() +
  labs(x = "1 - Specificity", y = "Sensitivity", color = "Model") +
  annotate("segment", x = 0, y = 0, xend = 1, yend = 1, linetype = 2, color = "gray") +
  ggtitle("ROC-Internal Validation")

print(g)
#ROC-test
g.list_xgbt<-ggroc(roc.list_xgbt,alpha=1,size=1,legacy.axes=TRUE,color="#FF4A46")
g.list_knnt<-ggroc(roc.list_knnt,alpha=1,size=1,legacy.axes=TRUE,color="#00BDcd")
g.list_gbmt<-ggroc(roc.list_gbmt,alpha=1,size=1,legacy.axes=TRUE,color="#008941")
g.list_enett<-ggroc(roc.list_enett,alpha=1,size=1,legacy.axes=TRUE,color="#993399")
g.list_logt<-ggroc(roc.list_logt,alpha=1,size=1,legacy.axes=TRUE,color="gold")
roc_list <- list(
  Enet = roc.list_enett,
  lightGBM = roc.list_gbmt,
  KNN = roc.list_knnt,
  XGBoost = roc.list_xgbt,
  LR=roc.list_logt
)

colors <- c(XGBoost = "#FF4A46", KNN = "#00BDcd", lightGBM = "#008941", Enet = "#993399",LR="gold")

roc_data <- do.call(rbind, lapply(names(roc_list), function(name) {
  roc_obj <- roc_list[[name]]
  data.frame(
    FPR = 1 - roc_obj$specificities, 
    TPR = roc_obj$sensitivities, 
    Model = name
  )
}))

g <- ggplot(roc_data, aes(x = FPR, y = TPR, color = Model)) +
  geom_line(size = 1, alpha = 1) +
  scale_color_manual(values = colors) + 
  theme_minimal() +
  labs(x = "1 - Specificity", y = "Sensitivity", color = "Model") +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), linetype = 2, color = "gray") +
  ggtitle("ROC-External Validation")

print(g)

#DCA-internal
library(dcurves)
#XGBoost
t1<-dat.train
t2<-dat.test
t1$XGBoost<-train_predictions
t2$XGBoost<-test_predictions
#KNN
t1$KNN<-knn_train
t2$KNN<-knn_test
#lightGBM
t1$lightGBM<-yt_pred_probs
t2$lightGBM<-y_pred_probs
#Enet
t1$Enet<-train_pred_probs
t2$Enet<-test_pred_probs
#LR
t1$LR<-gd
t2$LR<-gd2

dca_result <- dca(diagnosis ~ Enet + lightGBM + KNN + XGBoost+LR,
                  data = t1,
                  as_probability = c("Enet", "lightGBM", "KNN", "XGBoost","LR"))

p <- plot(dca_result) +
  labs(title = "DCA Curve - Internal Validation",
       x = "Threshold Probability",
       y = "Net Benefit") +
  theme_minimal() +
  theme(legend.position = "top") +
  geom_abline(intercept = mean(t1$diagnosis), slope = -1, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")


p <- p + scale_color_manual(values = c(XGBoost = "#FF4A46", KNN = "#00BDcd", lightGBM = "#008941", Enet = "#993399",LR="gold"))


p <- p + geom_line(size = 1) 

p <- p + scale_color_manual(values = c(XGBoost = "#FF4A46", KNN = "#00BDcd", lightGBM = "#008941", Enet = "#993399",LR="gold",
                                       treatAll = "black", treatNone = "black"),
                            labels = c(XGBoost = "XGBoost", KNN = "KNN", lightGBM = "lightGBM", Enet = "Enet",
                                       treatAll = "Treat All", treatNone = "Treat None"))

print(p)


dca_result <- dca(diagnosis ~ Enet + lightGBM + KNN + XGBoost+LR,
                  data = t2,
                  as_probability = c("Enet", "lightGBM", "KNN", "XGBoost","LR"))

p <- plot(dca_result) +
  labs(title = "DCA Curve - External Validation",
       x = "Threshold Probability",
       y = "Net Benefit") +
  theme_minimal() +
  theme(legend.position = "top") +
  geom_abline(intercept = mean(t2$diagnosis), slope = -1, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")


p <- p + scale_color_manual(values = c(XGBoost = "#FF4A46", KNN = "#00BDcd", lightGBM = "#008941", Enet = "#993399",LR="gold"))


p <- p + geom_line(size = 1)

p <- p + scale_color_manual(values = c(XGBoost = "#FF4A46", KNN = "#00BDcd", lightGBM = "#008941", Enet = "#993399",LR="gold",
                                       treatAll = "black", treatNone = "black"),
                            labels = c(XGBoost = "XGBoost", KNN = "KNN", lightGBM = "lightGBM", Enet = "Enet",LR="LR",
                                       treatAll = "Treat All", treatNone = "Treat None"))


print(p)
