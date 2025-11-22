
#install.packages("pROC")

library(pROC)
library(ggvenn)
setwd=("C:\\Users\\anxin\\Desktop\\integrated_model")
expFile="diffGeneExp.txt" 

lasso_genes <- scan("C:\\Users\\anxin\\Desktop\\integrated_model\\LASSO.gene.txt", what = character(), sep = "\n", quiet = TRUE)

rf_genes <- scan("C:\\Users\\anxin\\Desktop\\integrated_model\\rfGenes.txt", what = character(), sep = "\n", quiet = TRUE)


x<-list('LASSO genes'=lasso_genes,'RF genes'=rf_genes)
ggvenn(x,
       c("LASSO genes","RF genes"),
       show_percentage=FALSE,
       stroke_color="white",
       fill_color=c("#ffb2b2","#b2d4ec"),
       set_name_color=c("#ff0000","#1d6295"),
       text_size = 12)

geneRT=intersect(lasso_genes, rf_genes)
geneRT=as.data.frame(geneRT)


rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
y=gsub("(.*)\\_(.*)", "\\2", colnames(rt))
y=ifelse(y=="con", 0, 1)

for(x in as.vector(geneRT[,1])){
  roc1=roc(y, as.numeric(rt[x,]))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  
  library(ggplot2)
  g.list1 <- ggroc(roc1, alpha=1, size=1, legacy.axes=TRUE, color="#1B90CF")
  auc_value <- sprintf("%.4f", auc(roc1))
  plot <- g.list1 + 
    theme_update() +
    annotate(geom="segment", x=0, y=0, xend=1, yend=1, linetype=2) +
    annotate("text", x=0.7, y=0.20,
             label=paste0("AUC=", auc_value), colour="#5d6174", size=8) +
    ggtitle(paste0("ROC-", x))
  
  pdf(file=paste0("ROC.",x,".pdf"), width=5, height=5)
  print(plot)
  dev.off()
}

