

#install.packages("glmnet")


set.seed(123)
library(glmnet) 
library(pROC)
library(ggplot2)

rt=read.table("C:\\Users\\anxin\\Desktop\\integrated_model\\diffGeneExp.txt", header=T, sep="\t", check.names=F, row.names=1)
rt=t(rt)


x=as.matrix(rt)
y=gsub("(.*)\\_(.*)", "\\2", row.names(rt))
fit=glmnet(x, y, family = "binomial", alpha=1)
plot(fit,xvar="lambda",label=TRUE)
cvfit=cv.glmnet(x, y, family="binomial", alpha=1,type.measure='deviance',nfolds = 10)
pdf(file="cvfit.pdf",width=6,height=5.5)
plot(cvfit)
dev.off()



coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
lassoGene=row.names(coef)[index]
lassoGene=lassoGene[-1]

write.table(lassoGene, file="LASSO.gene.txt", sep="\t", quote=F, row.names=F, col.names=F)

three_genes <- c("PIRT", "PPP1R1C", "KRT83")


# Bootstrap

n_bootstrap <- 100
sample_size <- 0.8
gene_counts <- setNames(rep(0, length(three_genes)), three_genes)

boot_roc_data <- data.frame()


for (i in 1:n_bootstrap) {
  boot_indices <- sample(1:nrow(x), size = nrow(x)*sample_size, replace = TRUE)
  x_boot <- x[boot_indices, ]
  y_boot <- y[boot_indices]
  
  cvfit_boot <- cv.glmnet(x_boot, y_boot, family = "binomial", alpha = 1, 
                          type.measure = 'deviance', nfolds = 10)
  coef_boot <- coef(cvfit_boot, s = cvfit_boot$lambda.min)
  selected_genes <- row.names(coef_boot)[which(coef_boot != 0)]
  selected_genes <- selected_genes[selected_genes != "(Intercept)"]
  
  for (gene in three_genes) {
    if (gene %in% selected_genes) {
      gene_counts[gene] <- gene_counts[gene] + 1
    }
  }
  
  if (length(selected_genes) > 0) {  
    pred_prob <- predict(cvfit_boot, newx = x, s = cvfit_boot$lambda.min, type = "response")
    roc_obj <- roc(y, as.vector(pred_prob))
    temp_roc <- data.frame(
      specificity = roc_obj$specificities,
      sensitivity = roc_obj$sensitivities,
      bootstrap = i 
    )
    boot_roc_data <- rbind(boot_roc_data, temp_roc)
  }
  
  if (i %% 10 == 0) {
    cat(paste0(i, "resampling\n"))
  }
}


#ROC）

cvfit_full <- cv.glmnet(x, y, family = "binomial", alpha = 1, 
                        type.measure = 'deviance', nfolds = 10)
pred_prob_full <- predict(cvfit_full, newx = x, s = cvfit_full$lambda.min, type = "response")
roc_full <- roc(y, as.vector(pred_prob_full))
full_roc_data <- data.frame(
  specificity = roc_full$specificities,
  sensitivity = roc_full$sensitivities
)
full_auc <- round(auc(roc_full), 3)



pdf("Bootstrap_ROC_lasso.pdf", width=7, height=6)
ggplot() +
  geom_line(data = boot_roc_data, 
            aes(x = 1 - specificity, y = sensitivity, group = bootstrap), 
            color = "#a7b9d7", alpha = 0.3) +
  geom_line(data = full_roc_data, 
            aes(x = 1 - specificity, y = sensitivity), 
            color = "#576fa0", linewidth = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", alpha = 0.5) +
  annotate("text", x = 0.7, y = 0.2, 
           label = paste0("Overall AUC = ", full_auc), 
           color = "#576fa0", fontface = "bold") +
  labs(x = "1 - Specificity", 
       y = "Sensitivity", 
       title = "Bootstrap-ROC for LASSO") +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 12))
dev.off()

