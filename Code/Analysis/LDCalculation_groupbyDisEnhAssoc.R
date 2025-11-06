# install.packages("LDlinkR")
# R client for the LDlink web API (from NCI/NIH).
# Works directly with dbSNP rsIDs (e.g., rs4103200) and retrieves LD information using 1000 Genomes Project reference panels.
# Very convenient if you already have SNP IDs and want population-specific LD values (r², D′).
# Go to https://ldlink.nih.gov/?tab=apiaccess to register a token
library(LDlinkR)

ld_matrix <- LDmatrix(snps = c("rs4103200","rs9891119","rs8075676","rs8076051"),
                      pop = "ALL",      # European population
                      r2d = "r2",       # choose "r2" or "d"
                      token = "a10ee3ae4750")



# library(snpStats)
# # Works with genotype data (not just SNP IDs).
# # Can calculate LD (r², D′) between SNPs in a SnpMatrix object.
# 
# # geno is a SnpMatrix with individuals in rows, SNPs in columns
# ld_result <- ld(geno[, c("rs4103200","rs9891119","rs8075676","rs8076051")],
#                 stats = c("D.prime","R.squared"))


# # Install if not already
# # if(!require(LDheatmap)) install.packages("LDheatmap")
# library(LDheatmap)
# 
# # Suppose your LDmatrix result is called ld_matrix
# # The first column usually contains SNP IDs, so remove it
# ld_mat <- as.matrix(ld_matrix[,-1])
# rownames(ld_mat) <- ld_matrix$RS_number
# colnames(ld_mat) <- ld_matrix$RS_number
# 
# # Simple LD heatmap
# LDheatmap(ld_mat, genetic.distances = NULL, color = colorRampPalette(c("white","red"))(20))
# 
# library(reshape2)
# library(ggplot2)
# 
# # Reshape matrix to long format
# ld_long <- melt(ld_mat, varnames = c("SNP1","SNP2"), value.name = "r2")
# 
# # Heatmap
# ggplot(ld_long, aes(x = SNP1, y = SNP2, fill = r2)) +
#   geom_tile() +
#   scale_fill_gradient(low="white", high="red") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust=1)) +
#   labs(title = "LD heatmap (r²)", fill="r²")


# Set working directory to where your files are stored (adjust path as needed)
setwd("~/Manuscripts/98RWDisEnhPlus")


# ------------------------------------------------------------
# Keep one SNP per LD cluster based on lowest p-value
# Inputs:
#   ld_mat  = symmetric LD matrix (r² values, row/col = SNP IDs)
#   pvals   = named numeric vector of p-values (names = SNP IDs)
#   r2_thresh = LD threshold (default 0.2)
# Output:
#   A list with 'keep' (selected SNPs) and 'removed'
# ------------------------------------------------------------
prune_by_pvalue <- function(ld_mat, pvals, r2_thresh = 0.8) {
  stopifnot(all(rownames(ld_mat) == colnames(ld_mat)))
  snps <- rownames(ld_mat)
  remaining <- snps
  keep <- character()
  removed <- character()
  
  while (length(remaining) > 0) {
    # among remaining, pick SNP with smallest p-value
    pick <- names(which.min(pvals[remaining]))
    keep <- c(keep, pick)
    
    # find SNPs in high LD with 'pick' (including itself)
    cluster <- remaining[ld_mat[pick, remaining] >= r2_thresh]
    
    # mark the others for removal
    removed <- c(removed, setdiff(cluster, pick))
    
    # remove this cluster from further consideration
    remaining <- setdiff(remaining, cluster)
  }
  
  list(keep = keep, removed = removed)
}


topk=10

file <- paste0("Results/RWDisEnhPlus_predict_top", topk, "_Evidence.csv_group.csv")

data <- read.csv(file, header = TRUE)

library(dplyr)

data_byAssoc <- data %>%
  group_by(trait, eid) %>%
  summarise(
    hgnc        = paste(unique(hgnc), collapse = ", "),
    rsid        = paste(unique(rsid), collapse = ", "),
    consequence = paste(unique(consequence), collapse = ", "),
    pmid        = paste(unique(pmid), collapse = ", "),
    .groups = "drop"
  ) %>%
  mutate(Assoc = paste(trait, eid, sep = "_")) %>%
  select(Assoc, everything())


rowid2LDMatrix <- list()
rsid_list <- list()
pvals_list <- list()

for(i in 1:nrow(data_byAssoc)){
  print(i)
  txt = as.character(data_byAssoc[i,"rsid"])
  # Split by comma
  parts <- strsplit(txt, ",\\s*")[[1]]
  
  # Extract rsID and p-value using regex
  rsid <- sub("^(rs[0-9]+)\\s*\\(.*$", "\\1", parts)
  pval <- as.numeric(sub("^rs[0-9]+\\s*\\((.*)\\)$", "\\1", parts))
  # Create named numeric vector
  pvals <- setNames(pval, rsid)
  
  # print(rsid)
  # print(pvals)
  
  if(length(unique(rsid))<2){
    print(rsid)
    next
  }
  ld_matrix <- LDmatrix(snps = rsid,
                        pop = "ALL",
                        r2d = "r2",       # choose "r2" or "d"
                        token = "a10ee3ae4750")
  rowid2LDMatrix[[i]] = ld_matrix
  
  rsid_list[[i]] = rsid
  pvals_list[[i]] = pvals
}

LDthres = 0.2
for(i in 1:length(rowid2LDMatrix)){
  if(is.null(rowid2LDMatrix[[i]])) next
  cat("ROW",i)
  ld_matrix = rowid2LDMatrix[[i]]
  # print(ld_matrix)
  
  ld_mat <- as.matrix(ld_matrix[,-1])
  rownames(ld_mat) = colnames(ld_mat)
  print(ld_mat)
  
  pvals = pvals_list[[i]]
  pvals <- pvals[colnames(ld_mat)]
  print(pvals)
  
  res <- prune_by_pvalue(ld_mat, pvals, r2_thresh=LDthres)
  print(paste("Keep:",toString(res$keep)))    # SNPs kept (lowest p-value per cluster)
  print(paste("Remove:", toString(res$removed))) # SNPs removed
}

# # Example LD matrix (from LDmatrix)
# ld_mat <- matrix(c(
#   1.0, 0.85, 0.20, 0.10,
#   0.85, 1.0, 0.25, 0.05,
#   0.20, 0.25, 1.0, 0.90,
#   0.10, 0.05, 0.90, 1.0
# ), nrow=4, byrow=TRUE)
# rownames(ld_mat) <- colnames(ld_mat) <- c("rs4103200","rs9891119","rs8075676","rs8076051")
# 
# # Example p-values
# pvals <- c(rs4103200=1e-4, rs9891119=5e-3, rs8075676=1e-6, rs8076051=2e-5)
# 
# # Run pruning
# res <- prune_by_pvalue(ld_mat, pvals, r2_thresh=0.8)
# 
# res$keep    # SNPs kept (lowest p-value per cluster)
# res$removed # SNPs removed

