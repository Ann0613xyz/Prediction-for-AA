

#install.packages("ggplot2")
#install.packages("ggrepel")


library(dplyr)
library(ggplot2)
library(ggrepel)


logFCfilter=1           
adj.P.Val.Filter=0.05      




rt = read.table("C:\\Users\\anxin\\Desktop\\integrated_model\\all.txt", header=T, sep="\t", check.names=F)
Sig=ifelse((rt$adj.P.Val<adj.P.Val.Filter) & (abs(rt$logFC)>logFCfilter), ifelse(rt$logFC>logFCfilter,"Up","Down"), "Not")
#Volcano plot in train set
rt = mutate(rt, Sig=Sig)
p = ggplot(rt, aes(logFC, -log10(adj.P.Val)))+
  geom_point(aes(col=Sig))+
  scale_color_manual(values=c("green", "black","red"))+
  labs(title = " ")+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))

p1 = p + geom_label_repel(data=filter(rt, ((rt$adj.P.Val < adj.P.Val.Filter) & (abs(rt$logFC) > logFCfilter))),
                          box.padding=0.1, point.padding=0.1, min.segment.length=0.05,
                          size=1.5, aes(label=id)) + theme_bw()

print(p1, width=7, height=6.1)



#Volcano plot in test set
rt = read.table("C:\\Users\\anxin\\Desktop\\integrated_model\\testall.txt", header=T, sep="\t", check.names=F)
Sig=ifelse((rt$adj.P.Val<adj.P.Val.Filter) & (abs(rt$logFC)>logFCfilter), ifelse(rt$logFC>logFCfilter,"Up","Down"), "Not")


rt = mutate(rt, Sig=Sig)
p = ggplot(rt, aes(logFC, -log10(adj.P.Val)))+
  geom_point(aes(col=Sig))+
  scale_color_manual(values=c("green", "black","red"))+
  labs(title = " ")+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))

p1 = p + geom_label_repel(data=filter(rt, ((rt$adj.P.Val < adj.P.Val.Filter) & (abs(rt$logFC) > logFCfilter))),
                          box.padding=0.1, point.padding=0.1, min.segment.length=0.05,
                          size=1.5, aes(label=id)) + theme_bw()

print(p1, width=7, height=6.1)
