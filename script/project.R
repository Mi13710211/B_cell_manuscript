

#input the dataset
setwd("D:/analysis/project")
raw_data<-read.csv("D:/analysis/project/data_GSE255720.csv", header = TRUE, sep = ",")

library(tidyr) 
library(tidyverse)
library(dplyr)

#data wrangling
#change column names and remove additional rows and columns

colnames(raw_data)<- raw_data[1,]
clean_data<- raw_data[!(rownames(raw_data) == "Sample_type"), ]
colnames(clean_data)[2]<-"gene_symbol"
clean_data<- clean_data[,-1]
clean_data<- clean_data[,-2]
str(clean_data)

#remove not available in gene_symbol column
clean_data[clean_data==""]<-NA
clean_data <- drop_na(clean_data, gene_symbol)

#replace the gene_symbol with ENSEMBL_ID
ensembl_data<- read.csv("D:/analysis/project/ensembl_data.csv",header = TRUE, sep = "," )
clean_data_mix <- clean_data %>%
  left_join(ensembl_data, by = "gene_symbol")
any(is.na(clean_data_mix))
clean_data_mix <- drop_na(clean_data_mix, ENSEMBL.ID)
clean_data_mix <- clean_data_mix %>% select(-gene_symbol)

#determine ENSEMBL_ID as the row names
rownames(clean_data_mix)<- clean_data_mix$ENSEMBL.ID
data<-clean_data_mix %>% select(-ENSEMBL.ID)

# converted the count data values to integers and removed genes with a total count lower than 10
data <- data %>% mutate(across(everything(), as.numeric)) %>% round(0)
count<- data[which(rowSums(data)>10),]

``````````````````````````````````````````````````````````````
#B cell-related genes visualized in a heatmap 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")
library(edgeR)
install.packages("limma")
library(limma)
install.packages("RColorBrewer")
library(RColorBrewer)
if (!require("ggplot2")) install.packages("ggplot2")
library(ggplot2)
install.packages("pheatmap")
library(pheatmap)

#x is my count data in which columns are samples(UC, non.IBD) and rows are genes 
x<- clean_data
logcounts <- log2(x + 1)
target_genes <- c( "CD19"  ,  "CD79B" ,"CD79A", "SDC1", "PTPRC") # select target genes
mypalette <- brewer.pal(11, "RdYlBu")
condition<- factor(c( "UC"   ,    "UC"  ,   "non.IBD"  , "non.IBD", "non.IBD" ,"UC"  ,   "UC"   ,  "non.IBD", "UC"  ,   "non.IBD")) #determine condition for columns
morecols <- colorRampPalette(mypalette)
sorted_indices <- order(condition)
logcounts_sorted <- logcounts[target_genes, sorted_indices]
condition_sorted <- factor(condition[sorted_indices], levels = c("UC", "non.IBD"))
annotation_col <- data.frame(Group = condition_sorted)
rownames(annotation_col) <- colnames(logcounts_sorted)

pheatmap(logcounts_sorted,
         color = rev(morecols(50)),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,  
         main = "Expression of B cell-related genes",
         scale = "row",
         fontsize = 12,
         fontface_row = "bold",
         border_color = "darkgray",
         legend = TRUE,
         legend_breaks = seq(-2, 2, 0.5),
         annotation_legend = TRUE,
         legend_position = "left",
         annotation_legend_side = "left",
         gaps_col = NULL,
         fontsize_col = 12,
         cellwidth = 25,
         cellheight = 25,
         margins = c(10, 10),
         show_colnames = TRUE
)
``````````````````````````````````````````````````````````````
#Gene Set Enrichment Analysis (GSEA) visualized in a bar chart 
#Differential gene expression
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("DESeq2",force = TRUE )
library(DESeq2)
count<- data[which(rowSums(data)>10),]
condition<- factor(c( "UC"   ,    "UC"  ,   "non.IBD"  , "non.IBD", "non.IBD" ,"UC"  ,   "UC"   ,  "non.IBD", "UC"  ,   "non.IBD")) #determine condition for columns
coldata<- data.frame(row.names = colnames(count), condition)
dds<- DESeqDataSetFromMatrix (countData = count, colData = coldata, design = ~condition)
dds<-DESeq(dds)
# Compare UC vs non.IBD
res <- results(dds, contrast = c("condition", "UC", "non.IBD"))
sigs<-na.omit(res)
sigs<-sigs[sigs$padj<0.05,]

#GSEA
BiocManager:: install("clusterProfiler")
library(clusterProfiler)
BiocManager::install("AnnotationDbi")
library(AnnotationDbi)
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
res1<- res[res$baseMean>10,]
res1<-res1[order(-res1$log2FoldChange),]
gene_list<- res1$log2FoldChange
names(gene_list) <- rownames(res1)               
gene_list <- sort(gene_list, decreasing = TRUE)  
head(gene_list)

library(clusterProfiler)
library(org.Hs.eg.db)
library(grid) 

gse<-gseGO(gene_list,
           ont = "BP",
           keyType = "ENSEMBL",
           OrgDb = "org.Hs.eg.db", 
           eps = 1e-300)
#GSEA bar chart
target_description <- c("antimicrobial humoral response",
                        "immunoglobulin production",
                        "humoral immune response", 
                        "B cell mediated immunity",
                        "T cell proliferation",
                        "T cell differentiation")

# Filter the rows based on the target descriptions
filtered_data <- gse_data[gse_data$Description %in% target_description, ]

# Create the bar chart using the filtered data
library(ggplot2)
##LINE bar chart
ggplot(data = filtered_data, aes(y = reorder(Description, enrichmentScore))) +
  geom_segment(aes(x = 0, xend = enrichmentScore, yend = Description), 
               color = "red", size = 1.5) +  
  geom_point(aes(x = enrichmentScore), 
             color = "red", size = 6) +  
  geom_text(aes(x = max(enrichmentScore) + 0.3, 
                label = format(pvalue, scientific = TRUE)), 
            hjust = 0, size = 5, color = "black") +  
  annotate("text", 
           x = max(filtered_data$enrichmentScore) + 0.4, 
           y = nrow(filtered_data) + 0.5, 
           label = "P-value", 
           hjust = 0, size = 4.5, fontface = "bold", color = "black") +  
  labs(
    title = " GSEA for UC vs non-IBD patients",
    x = "Normalised Enrichment Score",
    y = "Gene Set Description"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 13,face = "bold", color = "black"), 
    axis.text.x = element_text(size = 12, color = "black"),  
    axis.title = element_text(size = 14, face = "bold"),     
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  scale_x_continuous(expand = c(0.1, 0.1)) +  
  coord_cartesian(clip = "off")

``````````````````````````````````````````````````````````````
# Immunoglobulin gene expression visualised in a Volcano plot 
install.packages("ggplot2")
library(ggplot2)
install.packages("ggrepel")
library(ggrepel)
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
install.packages("textshaping")
library(textshaping)
library(org.Hs.eg.db)
BiocManager::install("AnnotationDbi")
library(AnnotationDbi)
sigs<-na.omit(res)
sigs.df<- as.data.frame(sigs)
sigs.df
sigs.df$symbol<- mapIds(org.Hs.eg.db, keys = row.names(sigs.df), keytype = "ENSEMBL", column = "SYMBOL")
##indicate target genes:
target_genes <- c("IGHG1", "IGHG2", "IGHG3", "IGHG4","IGHA1", "IGHA2", "IGHM", "IGHD", "IGHE" )
sigs.df$label <- ifelse(sigs.df$symbol %in% target_genes, sigs.df$symbol, NA)

EnhancedVolcano(
  sigs.df,
  x = "log2FoldChange",
  y = "pvalue",
  lab =NA,  
  pCutoff = 1e-4,       
  FCcutoff = 2,         
  title = "Volcano Plot",
  subtitle = "Immunoglobulin Differential Genes Expression",
  caption = "Condition: UC vs non-IBD",
  pointSize = 4.0,      
  labSize = 6.0        
) +
  # Adjust label positions for target genes
  geom_text_repel(
    data = subset(sigs.df, symbol %in% target_genes),
    aes(x = log2FoldChange, y = -log10(pvalue), label = symbol),
    size=6, 
    box.padding = 0.5,       
    point.padding = 0.3,     
    segment.color = "yellow",
    segment.size = 0.5,      
    max.overlaps = Inf       
  )+
  # Add halo around target genes
  geom_point(
    data = subset(sigs.df, symbol %in% target_genes),
    aes(x = log2FoldChange, y = -log10(pvalue)),
    size = 1.5,             
    shape = 21,             
    fill = "transparent",    
    color = "black",         
    stroke = 2               
  ) 


``````````````````````````````````````````````````````````````

