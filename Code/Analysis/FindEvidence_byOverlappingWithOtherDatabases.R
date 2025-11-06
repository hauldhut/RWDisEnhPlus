#Step 1: Load Required Libraries and Data
# Install packages if not already installed
# install.packages(c("ggplot2", "dplyr", "tidyr", "ggrepel", "pheatmap"))
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
library(pheatmap)



# Set working directory to where your files are stored (adjust path as needed)
setwd("~/Manuscripts/98RWDisEnhPlus")

topk=10

topk_file <- paste0("Results/RWDisEnhPlus_predict_top", topk, ".txt")
evid_file <- paste0("Results/RWDisEnhPlus_predict_top", topk, "_Evidence.csv")
enh2EntrezID_file <- paste0("Data/enh2EntrezID.txt")

library(hash)


# Read the file topk_file
df_topk <- read.table(topk_file, sep = "\t", header = FALSE,
                      stringsAsFactors = FALSE, quote = "")

# Create a new hash
DOID2enh <- hash()
DOID2Name <- hash()

# Fill the hash
for (i in seq_len(nrow(df_topk))) {
  doid <- df_topk[i,1]
  key <- df_topk[i, 2]                       # disease name
  val <- as.character(df_topk[i, 3:ncol(df_topk)]) # enhancer regions (vector)
  val <- val[val != ""]                 # remove empty columns if any
  DOID2Name[[doid]] <- key
  DOID2enh[[doid]] <- val
}

length(DOID2enh)
length(DOID2Name)

enh2DiseaseDBDE_file <- paste0("Data/EDRelation_DiseaseEnhancer.csv")

data <- read.csv(enh2DiseaseDBDE_file, stringsAsFactors = FALSE)

# Initialize an empty hash
Disease2EnhDE <- hash()

# Loop through each row and append enhancers to corresponding disease key
for (i in seq_len(nrow(data))) {
  dis <- data$disease[i]
  enh <- data$enhancer[i]

  # If key exists, append; otherwise, create new vector
  if (has.key(dis, Disease2EnhDE)) {
    Disease2EnhDE[[dis]] <- c(Disease2EnhDE[[dis]], enh)
  } else {
    Disease2EnhDE[[dis]] <- c(enh)
  }
}

print(Disease2EnhDE)



# enh2DiseaseDB_file <- paste0("Data/EDRelation_ENdb.csv")
enh2DiseaseDB_file <- paste0("Data/EDRelation_EnDisease.csv")

data <- read.csv(enh2DiseaseDB_file, stringsAsFactors = FALSE)

# Initialize an empty hash
Disease2Enh <- hash()

# Loop through each row and append enhancers to corresponding disease key
for (i in seq_len(nrow(data))) {
  dis <- data$disease[i]
  enh <- data$enhancer[i]
  
  # If key exists, append; otherwise, create new vector
  if (has.key(dis, Disease2Enh)) {
    Disease2Enh[[dis]] <- c(Disease2Enh[[dis]], enh)
  } else {
    Disease2Enh[[dis]] <- c(enh)
  }
}

# Optional: print example entries
print(Disease2Enh)
Disease2Enh[["DOID:0014667"]]

# commonDisease = intersect(hash::keys(Disease2Enh), hash::keys(DOID2enh))
commonDisease = intersect(hash::keys(Disease2Enh), hash::keys(Disease2EnhDE))
commonDisease

# for(doid in commonDisease){
#   predEnh = DOID2enh[[doid]]
#   dbEnh = Disease2Enh[[doid]]
#   cat(doid, ":",length(predEnh),length(dbEnh),"\n")
#   commonEnh = intersect(predEnh, dbEnh)
#   if(length(commonEnh)!=0){
#     cat(doid, ":", commonEnh, "\n")
#   }
# }


calc_overlap_ratio <- function(region1, region2, method = c("union", "shorter", "longer")) {
  method <- match.arg(method)
  
  # Parse regions
  parse_region <- function(r) {
    parts <- strsplit(r, "[:-]")[[1]]
    list(chr = parts[1], start = as.numeric(parts[2]), end = as.numeric(parts[3]))
  }
  
  r1 <- parse_region(region1)
  r2 <- parse_region(region2)
  
  # Check chromosome match
  if (r1$chr != r2$chr) return(0)
  
  # Calculate overlap
  overlap_len <- max(0, min(r1$end, r2$end) - max(r1$start, r2$start))
  len1 <- r1$end - r1$start
  len2 <- r2$end - r2$start
  
  # Compute ratio
  if (overlap_len == 0) return(0)
  
  ratio <- switch(
    method,
    union   = overlap_len / (max(r1$end, r2$end) - min(r1$start, r2$start)),
    shorter = overlap_len / min(len1, len2),
    longer  = overlap_len / max(len1, len2)
  )
  return(ratio)
}

# Example:
calc_overlap_ratio("chr8:129543949-129554294", "chr8:129167819-129168726")

for(doid in commonDisease){
  predEnh = DOID2enh[[doid]]
  dbEnh = Disease2Enh[[doid]]
  
  for(pe in predEnh){
    for(de in dbEnh){
      ratio = calc_overlap_ratio(pe, de)
      # cat(doid, ":", pe,":", de,":",ratio, "\n")
      
      if(ratio>0){
        cat(doid, ":", pe,"and", de,":", ratio, "\n")
      }
    }
  }
}



