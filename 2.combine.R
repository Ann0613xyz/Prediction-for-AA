
library(limma)
library(sva)


rt1=read.csv("C:\\Users\\anxin\\Desktop\\integration\\GSE45512\\GSE45512annotated_probe_matrix.csv", header=T, check.names=F)
rt1=as.matrix(rt1)
rownames(rt1)=rt1[,1]
exp1=rt1[,2:ncol(rt1)]
dimnames=list(rownames(exp1),colnames(exp1))
data1=matrix(as.numeric(as.matrix(exp1)),nrow=nrow(exp1),dimnames=dimnames)
data1=avereps(data1)
data1=data1[rowMeans(data1)>0,]
data1=as.data.frame(data1)

rt2=read.csv("C:\\Users\\anxin\\Desktop\\integration\\GSE58573\\GSE58573annotated_probe_matrix.csv", header=T, check.names=F)
rt2=as.matrix(rt2)
rownames(rt2)=rt2[,1]
exp2=rt2[,2:ncol(rt2)]
dimnames=list(rownames(exp2),colnames(exp2))
data2=matrix(as.numeric(as.matrix(exp2)),nrow=nrow(exp2),dimnames=dimnames)
data2=avereps(data2)
data2=data2[rowMeans(data2)>0,]
data2=as.data.frame(data2)


rt3=read.csv("C:\\Users\\anxin\\Desktop\\integration\\GSE68801\\GSE68801annotated_probe_matrix.csv", header=T, check.names=F)
rt3=as.matrix(rt3)
rownames(rt3)=rt3[,1]
exp3=rt3[,2:ncol(rt3)]
dimnames=list(rownames(exp3),colnames(exp3))
data3=matrix(as.numeric(as.matrix(exp3)),nrow=nrow(exp3),dimnames=dimnames)
data3=avereps(data3)
data3=data3[rowMeans(data3)>0,]
data3=as.data.frame(data3)


rt4=read.csv("C:\\Users\\anxin\\Desktop\\integration\\GSE74761\\GSE74761annotated_probe_matrix.csv", header=T, check.names=F)
rt4=as.matrix(rt4)
rownames(rt4)=rt4[,1]
exp4=rt4[,2:ncol(rt4)]
dimnames=list(rownames(exp4),colnames(exp4))
data4=matrix(as.numeric(as.matrix(exp4)),nrow=nrow(exp4),dimnames=dimnames)
data4=avereps(data4)
data4=data4[rowMeans(data4)>0,]
data4=as.data.frame(data4)

rt5=read.csv("C:\\Users\\anxin\\Desktop\\integration\\GSE80342annotated_probe_matrix.csv", header=T, check.names=F)
rt5=as.matrix(rt5)
rownames(rt5)=rt5[,1]
exp5=rt5[,2:ncol(rt5)]
dimnames=list(rownames(exp5),colnames(exp5))
data5=matrix(as.numeric(as.matrix(exp5)),nrow=nrow(exp5),dimnames=dimnames)
data5=avereps(data5)
data5=data5[rowMeans(data5)>0,]
data5=as.data.frame(data5)

rt6=read.csv("C:\\Users\\anxin\\Desktop\\integration\\GSE148346\\GSE148346annotated_probe_matrix.csv", header=T, check.names=F)
rt6=as.matrix(rt6)
rownames(rt6)=rt6[,1]
exp6=rt6[,2:ncol(rt6)]
dimnames=list(rownames(exp6),colnames(exp6))
data6=matrix(as.numeric(as.matrix(exp6)),nrow=nrow(exp6),dimnames=dimnames)
data6=avereps(data6)
data6=data6[rowMeans(data6)>0,]
data6=as.data.frame(data6)

par(mar = c(7, 4, 4, 2))

boxplot(data1,outline=F,notch=T,las=2)
boxplot(data2,outline=F,notch=T,las=2)
boxplot(data3,outline=F,notch=T,las=2)
boxplot(data4,outline=F,notch=T,las=2)
boxplot(data5,outline=F,notch=T,las=2)
boxplot(data6,outline=F,notch=T,las=2)

data1=as.data.frame(normalizeBetweenArrays(data1))
data2=as.data.frame(normalizeBetweenArrays(data2))
data3=as.data.frame(normalizeBetweenArrays(data3))
data4=as.data.frame(normalizeBetweenArrays(data4))
data5=as.data.frame(normalizeBetweenArrays(data5))
data6=as.data.frame(normalizeBetweenArrays(data6))
boxplot(data1,outline=F,notch=T,las=2)
boxplot(data2,outline=F,notch=T,las=2)
boxplot(data3,outline=F,notch=T,las=2)
boxplot(data4,outline=F,notch=T,las=2)
boxplot(data5,outline=F,notch=T,las=2)
boxplot(data6,outline=F,notch=T,las=2)


library(dplyr)

#Gene count in data6 is different from other datasets
name1 <- intersect(rownames(data1),rownames(data2))
name2 <- intersect(name1,rownames(data3))
name3 <- intersect(name2,rownames(data4))
name4 <- intersect(name3,rownames(data5))
expr <- cbind(data1[name4,],data2[name4,],data3[name4,],data4[name4,],data5[name4,])
name5 <- intersect(name4, rownames(data6))

par(mar = c(7, 4, 4, 2))
boxplot(expr,outline=F,notch=T,las=2)

zero_genes <- apply(expr, 1, function(x) all(x == 0))
cat("Number of all-zero genes:", sum(zero_genes), "\n")


missing_values <- sum(is.na(expr))
cat("Number of missing values:", missing_values, "\n")

expr <- expr[!zero_genes, ]

expr[is.na(expr)] <- colMeans(expr, na.rm = TRUE)

batch <- c(rep(1, ncol(data1)), rep(2, ncol(data2)), rep(3, ncol(data3)), rep(4, ncol(data4)), rep(5, ncol(data5)))
combined_combat <- ComBat(expr, batch)
boxplot(combined_combat,outline=F,notch=T,las=2)

setwd="C:\\Users\\anxin\\Desktop\\integration"

write.csv(combined_combat,"combined_combat.csv")

# --------------------------
# PCA：data1 to data6
# --------------------------


#install.packages("missMDA")
library(missMDA)
library(FactoMineR)
library(factoextra)
clean_pca_data <- function(data) {
  data_mat <- as.matrix(data)
  data_mat[is.infinite(data_mat)] <- NA
  data_mat[is.na(data_mat)] <- colMeans(data_mat, na.rm = TRUE)
  zero_rows <- apply(data_mat, 1, function(x) all(x == 0))
  if(sum(zero_rows) > 0) {
    data_mat <- data_mat[!zero_rows, ]
  }
  return(as.data.frame(data_mat))
}

common_genes_all <- intersect(name4, rownames(data6)) 


pca_data1 <- clean_pca_data(data1[common_genes_all, ])
pca_data2 <- clean_pca_data(data2[common_genes_all, ])
pca_data3 <- clean_pca_data(data3[common_genes_all, ])
pca_data4 <- clean_pca_data(data4[common_genes_all, ])
pca_data5 <- clean_pca_data(data5[common_genes_all, ])
pca_data6 <- clean_pca_data(data6[common_genes_all, ])


all_data <- cbind(pca_data1, pca_data2, pca_data3, pca_data4, pca_data5, pca_data6)


all_data_t <- t(all_data)


remaining_na <- sum(is.na(all_data_t))
if(remaining_na > 0) {
  cat("detect", remaining_na, "missing values\n")
  ncp_est <- estim_ncpPCA(all_data_t, ncp.max = 5) 
  imputed_data <- imputePCA(all_data_t, ncp = ncp_est$ncp) 
  pca_result_all <- PCA(imputed_data$completeObs, scale.unit = FALSE, ncp = 5, graph = FALSE)
} else {
  pca_result_all <- PCA(all_data_t, scale.unit = FALSE, ncp = 5, graph = FALSE)
}


pca_df_all <- as.data.frame(pca_result_all$ind$coord[, 1:2])
pca_df_all$Dataset <- factor(
  c(rep("GSE45512", ncol(pca_data1)), 
    rep("GSE58573", ncol(pca_data2)),
    rep("GSE68801", ncol(pca_data3)),
    rep("GSE74761", ncol(pca_data4)),
    rep("GSE80342", ncol(pca_data5)),
    rep("GSE148346", ncol(pca_data6))),
  levels = c("GSE45512", "GSE58573", "GSE68801", "GSE74761", "GSE80342", "GSE148346")
)


eig_val_all <- get_eigenvalue(pca_result_all)
var1_all <- paste0("PC1 (", round(eig_val_all[1, 2], 1), "%)")
var2_all <- paste0("PC2 (", round(eig_val_all[2, 2], 1), "%)")



print(ggplot(pca_df_all, aes(x = Dim.1, y = Dim.2, color = Dataset)) +
        geom_point(size = 3, alpha = 0.7) +
        labs(title = "PCA of all datasets",
             x = var1_all, y = var2_all, color = "Dataset") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")))

# --------------------------
# PCA：train and test set
# --------------------------


common_genes_combat <- intersect(rownames(combined_combat), rownames(data6))


combat_pca_data <- clean_pca_data(combined_combat[common_genes_combat, ])
data6_pca_data <- clean_pca_data(data6[common_genes_combat, ])

combined_pca_data <- cbind(combat_pca_data, data6_pca_data)
combined_pca_t <- t(combined_pca_data)


pca_result_combined <- PCA(combined_pca_t, scale.unit = FALSE, ncp = 5, graph = FALSE)


pca_df_combined <- as.data.frame(pca_result_combined$ind$coord[, 1:2])
pca_df_combined$Group <- factor(
  c(rep("train set", ncol(combat_pca_data)), 
    rep("test set", ncol(data6_pca_data))),
  levels = c("train set", "test set")
)


eig_val_comb <- get_eigenvalue(pca_result_combined)
var1_comb <- paste0("PC1 (", round(eig_val_comb[1, 2], 1), "%)")
var2_comb <- paste0("PC2 (", round(eig_val_comb[2, 2], 1), "%)")



print(ggplot(pca_df_combined, aes(x = Dim.1, y = Dim.2, color = Group)) +
        geom_point(size = 3, alpha = 0.7) +
        labs(title = "PCA: train set vs test set",
             x = var1_comb, y = var2_comb, color = "Group") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")))


rt1 <- read.table("D:\\Deskstop\\integration\\allcontrol.txt", header = FALSE, sep = "\t", check.names = FALSE, fileEncoding = "UTF-8")
sampleName1 <- as.vector(rt1[, 1])
sampleName1 <- unique(sampleName1)


rt2 <- read.table("D:\\Deskstop\\integration\\allAA.txt", header = FALSE, sep = "\t", check.names = FALSE, fileEncoding = "UTF-8")
sampleName2 <- as.vector(rt2[, 1])
sampleName2 <- unique(sampleName2)

conData=combined_combat[,sampleName1]
AAData=combined_combat[,sampleName2]
combined_combat=cbind(conData,AAData)
conNum=ncol(conData)
AANum=ncol(AAData)

Type=c(rep("con",conNum),rep("AA",AANum))
design <- model.matrix(~0+factor(Type))
colnames(design) <- c("con","AA")
fit <- lmFit(combined_combat,design)
cont.matrix<-makeContrasts(con-AA,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)


logFCfilter=1             
adj.P.Val.Filter=0.05


allDiff=topTable(fit2,adjust='fdr',number=200000)
allDiffOut=rbind(id=colnames(allDiff),allDiff)
write.table(allDiffOut, file="all.txt", sep="\t", quote=F, col.names=F)

outData=rbind(id=paste0(colnames(combined_combat),"_",Type),combined_combat)
write.table(outData, file="normalize.txt", sep="\t", quote=F, col.names=F)

diffSig=allDiff[with(allDiff, (abs(logFC)>logFCfilter & adj.P.Val < adj.P.Val.Filter )), ]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut, file="diff.txt", sep="\t", quote=F, col.names=F)


diffGeneExp=combined_combat[row.names(diffSig),]
diffGeneExpOut=rbind(id=paste0(colnames(diffGeneExp),"_",Type), diffGeneExp)
write.table(diffGeneExpOut, file="diffGeneExp.txt", sep="\t", quote=F, col.names=F)


geneNum=50
diffSig=diffSig[order(as.numeric(as.vector(diffSig$logFC))),]
diffGeneName=as.vector(rownames(diffSig))
diffLength=length(diffGeneName)
hmGene=c()
if(diffLength>(2*geneNum)){
  hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
  hmGene=diffGeneName
}
hmExp=combined_combat[hmGene,]
Type=c(rep("Con",conNum),rep("AA",AANum))
names(Type)=colnames(combined_combat)
Type=as.data.frame(Type)
library(pheatmap)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 8,
         fontsize_row=4,
         fontsize_col=8)
