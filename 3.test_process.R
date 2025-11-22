
library(limma)
library(sva)

data6=read.csv("C:\\Users\\anxin\\Desktop\\integrated_model\\GSE148346\\GSE148346annotated_probe_matrix.csv", header=TRUE, check.names=F)
data6=as.matrix(data6)
rownames(data6)=data6[,1]
exp6=data6[,2:ncol(data6)]
dimnames=list(rownames(exp6),colnames(exp6))
data6=matrix(as.numeric(as.matrix(exp6)),nrow=nrow(exp6),dimnames=dimnames)
data6=avereps(data6)
data6=data6[rowMeans(data6)>0,]
data6=as.data.frame(data6)


setwd=("C:\\Users\\anxin\\Desktop\\integrated_model") 

rt1 <- read.table("C:\\Users\\anxin\\Desktop\\integrated_model\\testcontrol.txt", header = FALSE, sep = "\t", check.names = FALSE, fileEncoding = "UTF-8")
sampleName1 <- as.vector(rt1[, 1])
sampleName1 <- unique(sampleName1)


rt2 <- read.table("C:\\Users\\anxin\\Desktop\\integrated_model\\testAA.txt", header = FALSE, sep = "\t", check.names = FALSE, fileEncoding = "UTF-8")
sampleName2 <- as.vector(rt2[, 1])
sampleName2 <- unique(sampleName2)


conData=data6[,sampleName1]
AAData=data6[,sampleName2]
data6=cbind(conData,AAData)
conNum=ncol(conData)
AANum=ncol(AAData)

Type=c(rep("con",conNum),rep("AA",AANum))
design <- model.matrix(~0+factor(Type))
colnames(design) <- c("con","AA")
fit <- lmFit(data6,design)
cont.matrix<-makeContrasts(con-AA,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)


logFCfilter=1           
adj.P.Val.Filter=0.05


allDiff=topTable(fit2,adjust='fdr',number=200000)
allDiffOut=rbind(id=colnames(allDiff),allDiff)
write.table(allDiffOut, file="testall.txt", sep="\t", quote=F, col.names=F)


outData=rbind(id=paste0(colnames(data6),"_",Type),data6)
write.table(outData, file="testnormalize.txt", sep="\t", quote=F, col.names=F)

diffSig=allDiff[with(allDiff, (abs(logFC)>logFCfilter & adj.P.Val < adj.P.Val.Filter )), ]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut, file="testdiff.txt", sep="\t", quote=F, col.names=F)

 
diffGeneExp=data6[row.names(diffSig),]
diffGeneExpOut=rbind(id=paste0(colnames(diffGeneExp),"_",Type), diffGeneExp)
write.table(diffGeneExpOut, file="testdiffGeneExp.txt", sep="\t", quote=F, col.names=F)

#Check if the feature genes are present in the differentially expressed genes of the test set.

feature_genes <- c("KRT83","PPP1R1C","PIRT")

is_in_row_names <- feature_genes %in% rownames(diffGeneExpOut)
print(is_in_row_names)

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
hmExp=data6[hmGene,]
Type=c(rep("Con",conNum),rep("AA",AANum))
names(Type)=colnames(data6)
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
