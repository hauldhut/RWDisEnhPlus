#Step 1: Load Required Libraries and Data
# Install packages if not already installed
# install.packages(c("ggplot2", "dplyr", "tidyr", "ggrepel", "pheatmap"))
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
library(pheatmap)
library(biomaRt)


# Set working directory to where your files are stored (adjust path as needed)
setwd("~/Manuscripts/98RWDisEnhPlus/Results")

topk=10

file <- paste0("RWDisEnhPlus_predict_top", topk, "_Evidence.csv")

data <- read.csv(file, header = TRUE)

data <- data %>%
  separate(rsid, into = c("rsid", "p_value"), sep = " ", extra = "merge")


rsids <- unique(data$rsid)

print(length(rsids))

# Step 1: Map SNPs to Genes Using biomaRt
# Install Bioconductor packages
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install(c("clusterProfiler", "ReactomePA", "org.Hs.eg.db"))

predicted_snps <-rsids
  
ensembl <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")

# Query BioMart to map rsIDs to genes
snp_attributes <- c("refsnp_id", "ensembl_gene_stable_id", "associated_gene")
snp_to_gene <- getBM(attributes = snp_attributes,
                     filters = "snp_filter",
                     values = predicted_snps,
                     mart = ensembl)
saveRDS(snp_to_gene,paste0("snp_to_gene_top",topk,".rdata"))

###################
topk=10
snp_to_gene = readRDS(paste0("snp_to_gene_top",topk,".rdata"))
# View the mapping
head(snp_to_gene)

# Extract unique gene IDs (Ensembl IDs) or gene symbols
# Note: "associated_gene" may contain gene symbols, but we’ll use Ensembl IDs for pathway analysis
unique_genes <- unique(snp_to_gene$ensembl_gene_stable_id)
unique_genes <- unique_genes[unique_genes != ""]  # Remove empty entries

# If you prefer gene symbols, you can use "associated_gene" and split if multiple genes are listed
gene_symbols <- unique(unlist(strsplit(snp_to_gene$associated_gene, ",")))
gene_symbols <- gene_symbols[gene_symbols != ""]  # Remove empty entries

# Convert gene symbols to Entrez IDs (required for clusterProfiler)
library(org.Hs.eg.db)
gene_entrez <- mapIds(org.Hs.eg.db, 
                      keys = gene_symbols, 
                      column = "ENTREZID", 
                      keytype = "SYMBOL", 
                      multiVals = "first")
gene_entrez <- gene_entrez[!is.na(gene_entrez)]  # Remove unmapped genes

# Step 2: Perform Pathway Enrichment Analysis Using clusterProfiler
# The clusterProfiler package supports enrichment analysis for KEGG pathways, GO terms, and other databases. 

library(clusterProfiler)

# Perform KEGG pathway enrichment
kegg_enrich <- enrichKEGG(gene = gene_entrez,
                          organism = "hsa",  # Human (Homo sapiens)
                          pvalueCutoff = 0.05,  # Adjust p-value threshold as needed
                          qvalueCutoff = 0.2,   # Adjust q-value threshold as needed
                          minGSSize = 10,       # Minimum gene set size
                          maxGSSize = 500)      # Maximum gene set size

# View the results
head(as.data.frame(kegg_enrich))

# Step 3: Visualize and Interpret the Results
# You can visualize the enrichment results using plots like dot plots or bar plots. Here’s how to do it with clusterProfiler:
# Load ggplot2 for visualization
library(ggplot2)

kegg_sig <- kegg_enrich@result %>%
  filter(p.adjust <= 0.05)              
countSig = nrow(kegg_sig)

fontsize = 14
# Dot plot for KEGG enrichment
p_KEGG = dotplot(kegg_enrich, showCategory = countSig) + 
  # ggtitle("KEGG Pathway Enrichment for Involved SNPs") +
  theme(
    plot.margin = margin(t = 5, r = 0, b = 5, l = 0, unit = "pt"),  # Increase left margin further if needed
    axis.text.y = element_text(size = fontsize),  # Y-axis text size
    axis.title.x = element_text(size = fontsize),  # X-axis title font size (bold)
    axis.text.x = element_text(size = fontsize),   # X-axis tick labels font size
    legend.title = element_text(size = fontsize),  # Legend title font size (bold)
    legend.text = element_text(size = fontsize),    # Legend text font size
    # plot.title = element_text(size = fontsize + 2, face = "bold", hjust = 0.5),  # Increase
    aspect.ratio = 3  # Increase aspect ratio (height/width) to make the plot area narrower
  ) +
  scale_y_discrete(labels = function(x) x) +  # Ensure labels are not modified
  coord_cartesian(clip = "off")  # Prevent clipping of labels
print(p_KEGG)

ggsave(paste0("kegg_enrich@result_top",topk,".pdf"), p_KEGG, width = 12, height = 10)
ggsave(paste0("kegg_enrich@result_top",topk,".png"), p_KEGG, width = 12, height = 10)

