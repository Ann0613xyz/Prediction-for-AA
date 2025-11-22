library(biomaRt)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
probe_matrix <- read.table("C:\\Users\\anxin\\Desktop\\integrated_model\\GSE148346\\annotation\\GSE148346_series_matrix.txt", header = TRUE, row.names = 1)
platform_file <- read.table("C:\\Users\\anxin\\Desktop\\integrated_model\\GSE148346\\annotation\\GPL570-55999.txt", 
                            header = TRUE, sep = "\t", quote = "", 
                            stringsAsFactors = FALSE, comment.char = "#")
# Extract the correspondence between probe IDs and gene symbols
probe_gene_mapping <- platform_file[, c("ID", "Gene.Symbol")]
# Convert the row names of the probe matrix into a data frame and retain the content of the probe matrix.
probe_df <- data.frame(ID = rownames(probe_matrix), probe_matrix, stringsAsFactors = FALSE)

#Merge the data frames to replace probe names with gene symbols.
merged_data <- merge(probe_df, probe_gene_mapping, by = "ID", all.x = TRUE)

  

# Save the replaced probe matrix
write.csv(merged_data, "C:\\Users\\anxin\\Desktop\\integrated_model\\GSE148346\\annotation\\GSE148346annotated_probe_matrix.csv", row.names = TRUE)
