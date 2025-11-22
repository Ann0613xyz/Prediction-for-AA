

#install.packages('e1071')

#if (!requireNamespace("BiocManager", quietly = TRUE))
    #install.packages("BiocManager")
#BiocManager::install("preprocessCore")

setwd=("D:\\Deskstop\\integrated_model")  

inputData <- read.table("D:\\Deskstop\\integrated_model\\normalize.txt", header = TRUE, sep = "\t", check.names = TRUE)
LM22<-read.table("D:\\Deskstop\\integrated_model\\LM22.txt",header=T,sep='\t',row.names=1)


inputData <- inputData[!is.na(inputData[, 1]), ]

write.table(inputData, file = "temp_input_file.txt", sep = "\t", quote = FALSE, row.names = TRUE)
     
source("CIBERSORT.R") 
inputFile2 <- "temp_input_file.txt"

outTab=CIBERSORT(sig_matrix="LM22.txt", mixture_file=inputFile2, perm=1000, QN=TRUE)


outTab=outTab[outTab[,"P-value"]<0.05,]
outTab=as.matrix(outTab[,1:(ncol(outTab)-3)])
outTab=rbind(id=colnames(outTab),outTab)
write.table(outTab, file="CIBERSORT-Results.txt", sep="\t", quote=F, col.names=F)


