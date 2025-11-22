
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tibble)
library(tidyr)
library(broom)
library(limma)

data <- read.table("C:\\Users\\anxin\\Desktop\\integrated_model\\diffGeneExp.txt", header = TRUE, sep="\t")
data=na.omit(data)
rownames(data)<-data[,1]
data<-data[,-1]
genes <- c("KRT83", "PPP1R1C", "PIRT")
data <- data[genes, ]


data_long <- data %>%
  rownames_to_column(var = "gene") %>%
  tidyr::gather(key = "sample", value = "expression", -gene)


data_long$group <- ifelse(grepl("_con", data_long$sample), "Con", "AA")

expression_matrix <- as.matrix(data)

group <- ifelse(grepl("_con", colnames(data)), "Con", "AA")
design <- model.matrix(~0 + factor(group))
colnames(design) <- c("Con", "AA")

fit <- lmFit(expression_matrix, design)
cont.matrix <- makeContrasts(AA - Con, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

results <- topTable(fit2, adjust = 'fdr', number = Inf)
results_subset <- results[rownames(results) %in% genes, ]


results_subset$significaCone <- ifelse(results_subset$adj.P.Val < 0.005, "***",
                                      ifelse(results_subset$adj.P.Val < 0.01, "**",
                                             ifelse(results_subset$adj.P.Val < 0.05, "*", "")))

data_long$gene_group <- factor(paste(data_long$gene, data_long$group, sep = "-"), 
                               levels = c("KRT83-Con", "KRT83-AA", 
                                          "PPP1R1C-Con", "PPP1R1C-AA", 
                                          "PIRT-Con", "PIRT-AA"))

p <- ggplot(data_long, aes(x = gene_group, y = expression, fill = group,color=group)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white") +
  geom_jitter(width = 0.1, alpha = 0.5) +
  scale_fill_manual(values = c("Con" = "#6baed6", "AA" = "#fc9272")) +
  scale_color_manual(values = c("Con" = "#2171b5", "AA" = "#cb181d")) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  labs(x = "Feature Genes", y = "Expression", fill = "Group",color="Group") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

for (i in seq_along(genes)) {
  gene <- genes[i]
  sig <- results_subset[gene, "significaCone"]
  p <- p + annotate("text", x = i * 2 - 0.5, y = max(data_long$expression) * 1.1, 
                    label = sig, size = 5)
}

print(p)
