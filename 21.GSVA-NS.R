#BiocManager::install("GSVA")
library(GSVA)
library(limma)
library(ggplot2)
library(ggpubr)
library(clusterProfiler)
library(ggthemes)
library(ggprism)
setwd=("C:\\Users\\anxin\\Desktop\\integrated_model") 

expr_matrix <- read.table("C:\\Users\\anxin\\Desktop\\integrated_model\\normalize.txt",
                          sep = "\t")
row.names(expr_matrix)<-expr_matrix[,1]
expr_matrix<-expr_matrix[,-1]
colnames(expr_matrix)<-expr_matrix[1,]
expr_matrix<-expr_matrix[-1,]

group <- ifelse(grepl("_AA$", colnames(expr_matrix)), "AA", "Con")
expr_matrix<-as.matrix(expr_matrix)
class(expr_matrix) <- "numeric"

gmtFile="neural stress genesets.v2024.1.Hs.gmt"
sh <- read.gmt("C:\\Users\\anxin\\Desktop\\integrated_model\\neural stress genesets.v2024.1.Hs.gmt")
gene_sets<-sh
gene_sets<-split(gene_sets$gene, gene_sets$term)

param <- GSVA::gsvaParam(
  exprData = expr_matrix, 
  geneSets = gene_sets,
  kcdf = "Gaussian"       
)

gsva_scores <- GSVA::gsva(param, verbose = TRUE)
design <- model.matrix(~ 0 + factor(group, levels = c("AA", "Con")))
colnames(design) <- c("AA", "Con")
fit <- lmFit(gsva_scores, design)
contrast <- makeContrasts(AA_vs_Con = AA - Con, levels = design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
diff <- topTable(fit2, number = Inf, adjust.method = "BH")
dat_plot <- data.frame(id=row.names(diff), p=diff$P.Value, tvalue= diff$t)
dat_plot$group <- ifelse(dat_plot$tvalue >0 ,1,-1) 


dat_plot$g <- "Not"
dat_plot$g[ dat_plot$tvalue>0 & dat_plot$p < 0.05 ] <- "Up"
dat_plot$g[ dat_plot$tvalue<0 & dat_plot$p < 0.05 ] <- "Down"
table(dat_plot$g)
dat_plot$g <- factor(dat_plot$g, levels=c('Up','Down','Not'))


dat_plot$color <- ifelse(dat_plot$g=="Not", '#cccccc',"black")

dat_plot <- subset(dat_plot, g != "Not")


dat_plot$id <- factor(dat_plot$id)
dat_plot$g <- factor(dat_plot$g, levels=c('Up','Down'))

dat_plot <- dat_plot[order(dat_plot$tvalue,decreasing = T), ]
dat_plot$id <- factor(dat_plot$id,levels = rev(dat_plot$id))


dat_plot$lable_hjust <- ifelse(dat_plot$tvalue>0, 1, 0)


dat_plot$lable_xloc <- ifelse(dat_plot$tvalue>0, -0.05, 0.05)


t_up <- min(dat_plot[dat_plot$g=="Up","tvalue"]);t_up
t_down <- max(dat_plot[dat_plot$g=="Down","tvalue"]);t_down


p <- ggplot(data = dat_plot,aes(x = id, y = tvalue, fill = g)) +
  geom_col() +
  coord_flip() + 
  scale_fill_manual(values = c('Up'= '#fc9272','Down'='#6baed6')) +
  geom_hline(yintercept = c(t_down,t_up), color = 'white', linewidth = 0.5,lty='dashed') +
  xlab('') + 
  ylab('t value of GSVA score, AA versus Con') + 
  guides(fill="none") +
  theme_prism(border = T) +
  theme( axis.text.y = element_blank(),
         axis.ticks.y = element_blank() ) +
  geom_text(data = dat_plot, aes(x=id, y = lable_xloc, label = id, hjust = lable_hjust), 
            color=dat_plot$color, size = 2) 

p
setwd=("C:\\Users\\anxin\\Desktop\\integrated_model") 
write.csv(dat_plot,file="GSVA-neural stress.csv",quote=FALSE)
ggsave("GSVA_neural stress.pdf",p,device="pdf",width=11,height=10)

