library(limma)
library(randomForest)

set.seed(123)
setwd=("C:\\Users\\anxin\\Desktop\\integrated_model")
data=read.table("C:\\Users\\anxin\\Desktop\\integrated_model\\diffGeneExp.txt",head=TRUE,sep="\t",check.names=F)


first_col <- names(data)[1]


con_cols <- grep("_con$", names(data), value = TRUE)
aa_cols <- grep("_AA$", names(data), value = TRUE)

new_order <- c(first_col, con_cols, aa_cols)
data <- data[, new_order]

data <- data[!is.na(data[, 1]), ]


row.names(data) <- data[, 1]
data <- data[, -1]

genes=read.table("C:\\Users\\anxin\\Desktop\\integrated_model\\diff.txt",header=TRUE,sep="\t",check.names=F)

genes=genes[order(genes$adj.P.Val,decreasing=F),]

rownames(genes) <- genes[, 1] 


genes <- genes[, -1]

data=data[rownames(genes)[1:20],]

con_cols <- grep("_con$", names(data), value = TRUE)

aa_cols <- grep("_AA$", names(data), value = TRUE)

Type <- c(rep(1, length(con_cols)), rep(2, length(aa_cols)))

x=as.matrix(t(data))
y=Type

rf <- randomForest(x, as.factor(y), ntree = 500)

plot(rf,main="Random forest", lwd=2)

optionTrees=which.min(rf$err.rate[,1])
optionTrees
rf2=randomForest(x,as.factor(y),ntree=optionTrees)


importance=importance(x=rf2)

varImpPlot(rf2,main="")


rfGenes=importance[order(importance[,"MeanDecreaseGini"],decreasing=TRUE),]

rfGenes=names(rfGenes[rfGenes>2])

write.table(rfGenes,file="rfGenes.txt",sep="\t",quote=F,col.names=F,row.names=F)

#bootstrap
library(pROC)    
library(ggplot2) 
x_rf <- x 
y_rf <- as.factor(y)
three_genes <- c("PIRT", "PPP1R1C", "KRT83")
optionTrees <- optionTrees 



n_bootstrap <- 100 
sample_size <- 0.8 

gene_counts_rf <- setNames(rep(0, length(three_genes)), three_genes)

boot_roc_data_rf <- data.frame()


for (i in 1:n_bootstrap) {
  boot_indices <- sample(1:nrow(x_rf), size = nrow(x_rf)*sample_size, replace = TRUE)
  x_boot <- x_rf[boot_indices, ] 
  y_boot <- y_rf[boot_indices]
  
  rf_boot <- randomForest(x_boot, y_boot, ntree = optionTrees)
  
  imp_boot <- importance(rf_boot) 
  selected_genes_rf <- rownames(imp_boot)[imp_boot[, "MeanDecreaseGini"] > 2]
  
  for (gene in three_genes) {
    if (gene %in% selected_genes_rf) {
      gene_counts_rf[gene] <- gene_counts_rf[gene] + 1
    }
  }
  
  if (length(selected_genes_rf) > 0) { 
    pred_prob_boot <- predict(rf_boot, newdata = x_rf, type = "prob")[, 2]
    roc_obj_boot <- roc(y_rf, pred_prob_boot)
    temp_roc_rf <- data.frame(
      specificity = roc_obj_boot$specificities,
      sensitivity = roc_obj_boot$sensitivities,
      bootstrap = i
    )
    boot_roc_data_rf <- rbind(boot_roc_data_rf, temp_roc_rf)
  }

  if (i %% 10 == 0) {
    cat(paste0( i, "Resampling\n"))
  }
}


rf_full <- randomForest(x_rf, y_rf, ntree = optionTrees)

pred_prob_full_rf <- predict(rf_full, newdata = x_rf, type = "prob")[, 2]

roc_full_rf <- roc(y_rf, pred_prob_full_rf)
full_roc_data_rf <- data.frame(
  specificity = roc_full_rf$specificities,
  sensitivity = roc_full_rf$sensitivities
)
full_auc_rf <- round(auc(roc_full_rf), 3)



pdf("Bootstrap_ROC_RF.pdf", width = 7, height = 6)
ggplot() +
  geom_line(data = boot_roc_data_rf,
            aes(x = 1 - specificity, y = sensitivity, group = bootstrap),
            color = "#a7b9d7", alpha = 0.3) +

  geom_line(data = full_roc_data_rf,
            aes(x = 1 - specificity, y = sensitivity),
            color = "#576fa0", linewidth = 1.2) +

  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", alpha = 0.5) +

  annotate("text", x = 0.7, y = 0.2,
           label = paste0("Overall AUC = ", full_auc_rf),
           color = "#576fa0", fontface = "bold") +

  labs(x = "1 - Specificity",
       y = "Sensitivity",
       title = "Bootstrap-ROC for Random Forest") +

  xlim(0, 1) +
  ylim(0, 1) +

  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14), 
        axis.title = element_text(size = 12)) 
dev.off()

