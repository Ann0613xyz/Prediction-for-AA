

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")


library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

inputFile="testall.txt"      
gmtFile="c5.go.v7.4.symbols.gmt"  
setwd("C:\\Users\\anxin\\Desktop\\integrated_model")  

rt=read.table(inputFile, header=T, sep="\t", check.names=F)
logFC=as.vector(rt[,2])
names(logFC)=as.vector(rt[,1])
logFC=sort(logFC, decreasing=T)


gmt=read.gmt(gmtFile)


kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$p.adjust<0.05,]
write.table(kkTab,file="GSEA.test_GOresult.txt",sep="\t",quote=F,row.names = F)


termNum=10   
kkUp=kkTab[kkTab$NES>0,]
if(nrow(kkUp)>=termNum){
  kkUp$absES <- abs(kkUp$enrichmentScore)
  kkUp <- kkUp[order(kkUp$absES, decreasing = TRUE), ]
  showTerm=row.names(kkUp)[1:termNum]
  gseaplot <- gseaplot2(kk, showTerm, 
                        base_size = 8, 
                        title = "Enriched in AA",
                        rel_heights = c(2.5, 0.5, 1))
  pdf(file="GSEA.testGO_AA_Top10.pdf", width=10, height=8)
  print(gseaplot)
  dev.off()
}



termNum=10 
kkDown=kkTab[kkTab$NES<0,]
if(nrow(kkDown)>=termNum){
  kkDown$absES <- abs(kkDown$enrichmentScore)
  kkDown <- kkDown[order(kkDown$absES, decreasing = TRUE), ]
  showTerm=row.names(kkDown)[1:termNum]
  gseaplot=gseaplot2(kk, 
                     showTerm, 
                     base_size=8, 
                     title="Enriched in Con",
                     rel_heights = c(2.5, 0.5, 1))
  pdf(file="GSEA.testGO_con_Top10.pdf", width=10, height=8) 
  print(gseaplot)
  dev.off()
}
