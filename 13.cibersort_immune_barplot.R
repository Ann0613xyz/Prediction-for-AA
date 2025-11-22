
#install.packages("corrplot")


library(corrplot)
inputFile="CIBERSORT-Results.txt" 
setwd=("C:\\Users\\anxin\\Desktop\\integrated_model") 


rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)


if (any(is.na(rt))) {
  cat("Missing values existing.\n")
}

rt = rt[, apply(rt, 2, sd) != 0]


con=grepl("_con", rownames(rt), ignore.case=T)
AA=grepl("_AA", rownames(rt), ignore.case=T)
conData=rt[con,]
AAData=rt[AA,]
conNum=nrow(conData)
AANum=nrow(AAData)
data=t(rbind(conData,AAData))

data_long <- reshape2::melt(data)
data_long$group <- ifelse(grepl("_con", data_long$Var1, ignore.case = TRUE), "Con", "AA")
pdf(file="cibersort_barplot.pdf", width=16, height=9)
col <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", 
                          "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))(nrow(data))
par(mar=c(5, 6, 4, 15), xpd=TRUE)
a1 <- barplot(data, col=col, yaxt="n", ylab="Relative Percent", xaxt="n", cex.lab=1.8)
a2 <- axis(2, tick=FALSE, labels=FALSE)
axis(2, at=a2, labels=paste0(a2 * 100, "%"))
par(srt=0, xpd=TRUE)
rect(xleft = a1[1], ybottom = -0.01, xright = a1[conNum], ytop = -0.06, col="lightgreen", border="darkgreen")
text(a1[conNum] / 2, -0.035, "Con", cex=1.8, col="darkgreen")
rect(xleft = a1[conNum], ybottom = -0.01, xright = a1[length(a1)], ytop = -0.06, col="lightcoral", border="darkred")
text((a1[length(a1)] + a1[conNum]) / 2, -0.035, "AA", cex=1.8, col="darkred")
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),col=col,pch=15,bty="n",cex=1.2)
dev.off()

conData <- conData[, apply(conData, 2, sd) != 0]
AAData <- AAData[, apply(AAData, 2, sd) != 0]

conCor <- cor(conData, use = "complete.obs", method = "pearson")
conTestRes <- cor.mtest(conData, conf.level = 0.95)
par(mfrow=c(2,3))

corrplot(conCor,method='circle',
         tl.cex=0.8)  
corrplot(
  conCor,
  method='color',
  type='upper',
  add=TRUE,
  tl.pos='n',
  cl.pos='n',
  diag=F,
  p.mat=conTestRes$p,
  sig.level=c(0.001,0.01,0.05),
  pch.cex=1.5,
  insig='label_sig'
)  

AACor <- cor(AAData, use = "complete.obs", method = "pearson")
AATestRes <- cor.mtest(AAData, conf.level = 0.95)


corrplot(AACor,method='circle',
         tl.cex=0.8) 
corrplot(
  AACor,
  method='color',
  type='upper',
  add=TRUE,
  tl.pos='n',
  cl.pos='n',
  diag=F,
  p.mat=AATestRes$p,
  sig.level=c(0.001,0.01,0.05),
  pch.cex=1.5,
  insig='label_sig'
)  
