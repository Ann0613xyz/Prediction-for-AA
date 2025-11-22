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
library(dplyr)
library(tibble)

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

dat.train=train[,feature_genes]

dat.train <- cbind(rownames(dat.train), dat.train)

colnames(dat.train)[1] <- "diagnosis"
dat.train[, 1] <- sub(".*_", "", dat.train[, 1])
dat.train[, 1][dat.train[, 1] == "con"] <- 0
dat.train[, 1][dat.train[, 1] == "AA"] <- 1

immune_cell=read.table("C:\\Users\\anxin\\Desktop\\integrated_model\\CIBERSORT-Results.txt", header=TRUE, sep="\t", check.names=F)
columns_to_keep <- c("T cells CD4 memory resting", 
                     "Neutrophils", 
                     "Dendritic cells activated", 
                     "Macrophages M1", 
                     "T cells gamma delta")

immune_cell_final <- immune_cell %>%
  column_to_rownames(var = "id") %>%
  select(all_of(columns_to_keep))
immune_cell_final<-na.omit(immune_cell_final)

immune_cell_joined <- immune_cell_final %>%
  rownames_to_column(var = "Sample_ID")


dat_train_joined <- dat.train %>%
  rownames_to_column(var = "Sample_ID")

merged_data <- dat_train_joined %>%
  inner_join(immune_cell_joined, by = "Sample_ID")


merged_data_final <- merged_data %>%
  column_to_rownames(var = "Sample_ID")

merged_data_final <- merged_data_final %>%
  mutate(across(-diagnosis, as.numeric))

min_max_scale<-function(x){
  (x-min(x))/(max(x)-min(x))
}

immune.train=merged_data_final%>%
  mutate_if(.predicate=is.numeric,
            .funs=min_max_scale)%>%
  as.data.frame()

#XGboost-SHAP
library(data.table)
library(xgboost)
library(Matrix)
library(caTools)

x_train=immune.train[,-1]
x_train <- data.frame(lapply(x_train, as.numeric))
y_train=immune.train[,1]

immune.train <- na.omit(immune.train)
x_train<-na.omit(x_train)


immune.train[, 1] <- as.character(immune.train[, 1])
immune.train[, 1] <- as.numeric(immune.train[, 1])


immune.train <- na.omit(immune.train)
x_train <- immune.train[, -1] 
y_train <- immune.train[, 1] 

x_train <- as.matrix(data.frame(lapply(x_train, as.numeric)))
y_train <- as.numeric(as.character(y_train))


dtrain <- xgb.DMatrix(data = x_train, label = y_train)
params<-list(objective="binary:logistic",eval_metric="logloss",eta=0.1,max_depth=3)
nrounds=100

xgb_model_final <- xgboost(
  params = params,
  data = dtrain,
  nrounds = nrounds
)

train_predictions<-predict(xgb_model_final,newdata=dtrain)
train_predictions1<-ifelse(train_predictions>0.5,1,0)

#SHAP
library(shapviz)
for (i in 2:ncol(immune.train)) {
  immune.train[[i]] <- as.numeric(immune.train[[i]])
}

original_colnames <- colnames(immune.train)
feature_colnames <- original_colnames[original_colnames != "diagnosis"]
renamed_features <- gsub(" ", ".", feature_colnames)
new_colnames <- colnames(immune.train)
for (i in seq_along(new_colnames)) {
  if (new_colnames[i] %in% feature_colnames) {
    new_colnames[i] <- gsub(" ", ".", new_colnames[i])
  }
}
colnames(immune.train) <- new_colnames

shap<-shapviz(xgb_model_final,X_pred=data.matrix(immune.train[,-1]))

setwd="C:\\Users\\anxin\\Desktop\\integrated_model"
pdf(file="immunecell_genes_SHAP_interaction.pdf", width=16, height=16)
shp_i<-shapviz(xgb_model_final,
               X_pred=data.matrix(immune.train[,-1]),
               interactions=TRUE)
sv_interaction(shp_i)+
  theme(axis.text.x=element_text(angle=45,
                                 vjust=1,
                                 hjust=1))
dev.off()

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
          file = "immune_included_SHAP_Interaction_Quantification.csv", 
          row.names = FALSE)

