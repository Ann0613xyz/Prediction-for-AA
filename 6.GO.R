

# install.packages("colorspace")
# install.packages("stringi")
# install.packages("ggplot2")
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("DOSE")
# BiocManager::install("clusterProfiler")
# BiocManager::install("enrichplot")


library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

pvalueFilter=0.05      
qvalueFilter=0.05     


colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}

setwd=("C:\\Users\\anxin\\Desktop\\integrated_model")      
rt=read.table("C:\\Users\\anxin\\Desktop\\integrated_model\\diff.txt", header=T, sep="\t", check.names=F)  


genes=as.vector(rt[,1])
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]


kk=enrichGO(gene=gene, OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]

write.table(GO, file="GO.txt", sep="\t", quote=F, row.names = F)


showNum=15

kk_df <- as.data.frame(kk)
kk_sorted_df <- kk_df[order(kk_df$qvalue), ]

kk_filtered_df <- kk_sorted_df[(kk_sorted_df$pvalue < pvalueFilter & kk_sorted_df$qvalue < qvalueFilter), ]

if(nrow(kk_filtered_df) > showNum) {
  kk_filtered_df <- kk_filtered_df[1:showNum, ]
}

kk_filtered <- kk
kk_filtered@result <- kk_filtered_df


kk_filtered@result$Description <- gsub(" ", "\u00A0", kk_filtered@result$Description)

bar <- barplot(kk_filtered, drop = TRUE, showCategory = showNum, split = "ONTOLOGY", color = colorSel) +
  theme(axis.text.x = element_text(size = 10, hjust = 1),
        axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12), 
        plot.title = element_text(size = 14, face = "bold"),
        plot.margin = unit(c(1, 1, 1, 1), "lines"),
        strip.text = element_text(size = 8), 
        panel.spacing = unit(0.2, "lines"))+
  facet_grid(ONTOLOGY ~ ., scale = 'free', space = 'free_y') 


print(bar)
ggsave("GO_barplot.pdf", bar, width = 12, height = 8, dpi = 300) 


bub <- dotplot(kk_filtered, showCategory = showNum, orderBy = "GeneRatio", split = "ONTOLOGY", color = colorSel) +
  theme(axis.text.x = element_text(size = 10, hjust = 1), 
        axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"),
        plot.margin = unit(c(1, 1, 1, 1), "lines"),
        strip.text = element_text(size = 8),
        panel.spacing = unit(0.2, "lines"))+
  facet_grid(ONTOLOGY ~ ., scale = 'free', space = 'free_y')


print(bub)
ggsave("GO_bubplot.pdf", bub, width = 12, height = 8, dpi = 300) 
