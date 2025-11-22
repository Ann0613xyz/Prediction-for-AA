
#install.packages("vioplot")


library(vioplot) 
inputFile="CIBERSORT-Results.txt" 
setwd=("C:\\Users\\anxin\\Desktop\\integrated_model") 

rt=read.table(inputFile, header=TRUE, sep="\t", check.names=F, row.names=1)


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
rt=rbind(conData,AAData)


outTab <- data.frame()
pdf(file = "cibersort vioplot.pdf", height = 8, width = 13)
par(las = 1, mar = c(10, 6, 3, 3))

x <- 1:ncol(rt)
y <- 1:ncol(rt)

plot(x, y,
     xlim = c(0, 63), ylim = c(min(rt), max(rt) + 0.05),
     main = "", xlab = "", ylab = "Fraction",
     pch = 21,
     col = "white",
     xaxt = "n")

for (i in 1:ncol(rt)) {
  if (sd(rt[1:conNum, i]) == 0) {
    rt[1, i] <- 0.00001
  }
  if (sd(rt[(conNum + 1):(conNum + AANum), i]) == 0) {
    rt[(conNum + 1), i] <- 0.00001
  }
  
  conData <- rt[1:conNum, i]
  AAData <- rt[(conNum + 1):(conNum + AANum), i]
  
  vioplot(conData, at = 3*(i-1), lty = 1, add = TRUE, col = '#6baed6')
  vioplot(AAData, at = 3*(i-1)+1, lty = 1, add = TRUE, col = '#fc9272')
  
  wilcoxTest <- wilcox.test(conData, AAData)
  p <- wilcoxTest$p.value
  
  if (p < 0.05) {
    cellPvalue <- cbind(Cell = colnames(rt)[i], pvalue = p)
    outTab <- rbind(outTab, cellPvalue)
  }
  
  mx <- max(c(conData, AAData))
  lines(c(x = 3*(i-1) + 0.2, x = 3*(i-1) + 0.8), c(mx, mx))
  

  text(x = 3*(i-1) + 0.5, y = mx + 0.02, 
       labels = ifelse(p < 0.001, "p<0.001", paste0("p=", sprintf("%.03f", p))), 
       cex = 0.8)

  sig_symbol <- ifelse(p < 0.001, "***",
                       ifelse(p < 0.01, "**",
                              ifelse(p < 0.05, "*", "")))
  if (sig_symbol != "") {
    text(x = 3*(i-1) + 0.5, y = mx + 0.05, 
         labels = sig_symbol, 
         cex = 1.2, col = "black")
  }
}

legend("topright", 
       c("Con", "AA"),
       lwd = 3, bty = "n", cex = 1,
       col = c("#6baed6", "#fc9272"))

text(seq(1, 64, 3), -0.05, xpd = NA, labels = colnames(rt), cex = 1, srt = 45, pos = 2)

dev.off()

write.table(outTab,file="immuneDiff.xls",sep="\t",row.names=F,quote=F)
