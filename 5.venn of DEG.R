
#install.packages("pROC")

library(pROC)
library(ggvenn)
setwd=("C:\\Users\\anxin\\Desktop\\integrated_model")
train=read.table("C:\\Users\\anxin\\Desktop\\integrated_model\\diff.txt",sep = "\t")
test=read.table("C:\\Users\\anxin\\Desktop\\integrated_model\\testdiff.txt",sep = "\t")
new_colnames <- as.character(train[1, ])
train <- train[-1, ]
colnames(train) <- new_colnames
train$logFC <- as.numeric(as.character(train$logFC))
train$id <- as.character(train$id)
trainup <- train$id[train$logFC > 0]
traindown <- train$id[train$logFC < 0]

new_colnames <- as.character(test[1, ])
test <- test[-1, ]
colnames(test) <- new_colnames
test$logFC <- as.numeric(as.character(test$logFC))
test$id <- as.character(test$id)
testup <- test$id[test$logFC > 0]
testdown <- test$id[test$logFC < 0]

x<-list('up-regulated DEGs\n in train'=trainup,'up-regulated DEGs\n in test'=testup)
ggvenn(x,
       c("up-regulated DEGs\n in train","up-regulated DEGs\n in test"),
       show_percentage=FALSE,
       stroke_color="white",
       fill_color=c("#ffb2b2","#b2d4ec"),
       set_name_color=c("#ff0000","#1d6295"),
       text_size = 12,
       set_name_size = 5)
y<-list('down-regulated DEGs\n in train'=traindown,'down-regulated DEGs\n in test'=testdown)
ggvenn(y,
       c("down-regulated DEGs\n in train","down-regulated DEGs\n in test"),
       show_percentage=FALSE,
       stroke_color="white",
       fill_color=c("#ffb2b2","#b2d4ec"),
       set_name_color=c("#ff0000","#1d6295"),
       text_size = 12,
       set_name_size = 5)
