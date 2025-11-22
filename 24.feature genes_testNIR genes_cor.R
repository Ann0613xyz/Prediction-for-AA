setwd=("C:\\Users\\anxin\\Desktop\\integrated_model")

expr_matrix <- read.table("C:\\Users\\anxin\\Desktop\\integrated_model\\testnormalize.txt",
                          sep = "\t")
row.names(expr_matrix)<-expr_matrix[,1]
expr_matrix<-expr_matrix[,-1]
colnames(expr_matrix)<-expr_matrix[1,]
expr_matrix<-expr_matrix[-1,]
expr_matrix_numeric <- apply(expr_matrix, 2, as.numeric)
rownames(expr_matrix_numeric) <- rownames(expr_matrix)
colnames(expr_matrix_numeric) <- colnames(expr_matrix)
expr_matrix<-t(expr_matrix_numeric)
feature_genes <- c("KRT83", "PPP1R1C", "PIRT")


library(clusterProfiler)
gmtFile="GOBP_NEUROINFLAMMATORY_RESPONSE.v2024.1.Hs.gmt"
sh <- read.gmt(gmtFile)
gene_sets<-unique(sh[,2])
valid_genes <- c()

for (gene in gene_sets) {
  if (gene %in% colnames(expr_matrix)) {
    valid_genes <- c(valid_genes, gene)
  }
}

cor_matrix <- cor(expr_matrix[, c(feature_genes, valid_genes)])
subset_cor_matrix <- cor_matrix[feature_genes, valid_genes]
library(pheatmap)
breaks <- seq(-1, 1, length.out = 100)
pdf("feature genes testNIRcor_heatmap.pdf", width = 15, height = 3)
pheatmap(subset_cor_matrix,
         breaks=breaks,
         scale = "none",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row=10,
         fontsize_col=8,
         angle_col=45,
         annotation_legend = TRUE,
         legend_labels = list(col = "ρ"),
         color = colorRampPalette(c("blue", "white", "red"))(100))
dev.off()
