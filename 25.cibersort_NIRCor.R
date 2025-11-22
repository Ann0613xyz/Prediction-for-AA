

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("reshape2")
#install.packages("ggpubr")
#install.packages("ggExtra")


library(limma)
library(reshape2)
library(ggpubr)
library(ggExtra)
library(dplyr)
library(clusterProfiler)
expFile="normalize.txt"            
immFile="CIBERSORT-Results.txt"    
setwd=("C:\\Users\\anxin\\Desktop\\integrated_model")   

gmtFile="GOBP_NEUROINFLAMMATORY_RESPONSE.v2024.1.Hs.gmt"
sh <- read.gmt(gmtFile)
gene_sets<-unique(sh[,2])

#train
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt <- rt[!is.na(rt[, 1]), ]
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)


immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)


immune <- immune[complete.cases(immune), ]


immune <- immune[, apply(immune, 2, sd, na.rm=TRUE) != 0]
corr_matrix <- matrix(NA, nrow = length(gene_sets), ncol = ncol(immune))
rownames(corr_matrix) <- gene_sets
colnames(corr_matrix) <- colnames(immune)

for (i in seq_along(gene_sets)) {
  gene <- gene_sets[i]
  
  if (gene %in% rownames(data)) {

    current_data <- t(data[gene, , drop = F])
    current_data <- as.data.frame(current_data)
    
    samele <- intersect(row.names(immune), row.names(current_data))
    if (length(samele) >= 3) {
      rt <- cbind(immune[samele, ], current_data[samele, ])
      colnames(rt)[ncol(rt)] <- gene
      
      for (j in seq_along(colnames(immune))) {
        cell <- colnames(immune)[j]
        x <- as.numeric(rt[, gene])
        y <- as.numeric(rt[, cell])
        
        non_na_indices <- complete.cases(x, y)
        x <- x[non_na_indices]
        y <- y[non_na_indices]
        
        if (length(x) >= 3) {
          if (sd(y) == 0) {
            y[1] <- 0.00001
          }
          cor_result <- cor.test(x, y, method = "spearman")
          corr_matrix[i, j] <- cor_result$estimate
        } else {
          corr_matrix[i, j] <- NA
        }
      }
    } else {
      corr_matrix[i, ] <- NA
    }
  } else {
    corr_matrix[i, ] <- NA
  }
}
corr_matrix <- corr_matrix[!apply(is.na(corr_matrix), 1, any), ]
corr_matrix <- corr_matrix[, !apply(is.na(corr_matrix), 2, any)]
corr_matrix<-t(corr_matrix)
library(pheatmap)
breaks <- seq(-1, 1, length.out = 100)
pdf("immune_NIRcor_heatmap.pdf", width = 15, height = 4)
pheatmap(corr_matrix, 
         scale = "none", 
         color = colorRampPalette(c("blue", "white", "red"))(100), 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean",
         annotation_legend = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row=10,
         fontsize_col=8,
         angle_col=45,
         breaks=breaks)
dev.off()

write.csv(corr_matrix,file="NI immune cor.csv",quote=FALSE)


#test
expFile="testnormalize.txt"           
immFile="CIBERSORT-testResults.txt"  
setwd=("C:\\Users\\anxin\\Desktop\\integrated_model") 



rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt <- rt[!is.na(rt[, 1]), ]
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)


immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)


immune <- immune[complete.cases(immune), ]


immune <- immune[, apply(immune, 2, sd, na.rm=TRUE) != 0]
corr_matrix <- matrix(NA, nrow = length(gene_sets), ncol = ncol(immune))
rownames(corr_matrix) <- gene_sets
colnames(corr_matrix) <- colnames(immune)

for (i in seq_along(gene_sets)) {
  gene <- gene_sets[i]
  
  if (gene %in% rownames(data)) {
    current_data <- t(data[gene, , drop = F])
    current_data <- as.data.frame(current_data)
    
    samele <- intersect(row.names(immune), row.names(current_data))
    if (length(samele) >= 3) {
      rt <- cbind(immune[samele, ], current_data[samele, ])
      colnames(rt)[ncol(rt)] <- gene
      
      for (j in seq_along(colnames(immune))) {
        cell <- colnames(immune)[j]
        x <- as.numeric(rt[, gene])
        y <- as.numeric(rt[, cell])
        
        non_na_indices <- complete.cases(x, y)
        x <- x[non_na_indices]
        y <- y[non_na_indices]
        
        if (length(x) >= 3) {
          if (sd(y) == 0) {
            y[1] <- 0.00001
          }
          cor_result <- cor.test(x, y, method = "spearman")
          corr_matrix[i, j] <- cor_result$estimate
        } else {
          corr_matrix[i, j] <- NA
        }
      }
    } else {
      corr_matrix[i, ] <- NA
    }
  } else {
    corr_matrix[i, ] <- NA
  }
}
corr_matrix <- corr_matrix[!apply(is.na(corr_matrix), 1, any), ]
corr_matrix <- corr_matrix[, !apply(is.na(corr_matrix), 2, any)]

corr_matrix<-t(corr_matrix)
library(pheatmap)
breaks <- seq(-1, 1, length.out = 100)
pdf("immune_testNIRcor_heatmap.pdf", width = 15, height = 4)
pheatmap(corr_matrix, 
         scale = "none", 
         color = colorRampPalette(c("blue", "white", "red"))(100), 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean",
         annotation_legend = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row=10,
         fontsize_col=8,
         angle_col=45,
         breaks=breaks)
dev.off()
