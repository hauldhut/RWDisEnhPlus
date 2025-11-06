library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)   # for mapping Entrez -> Symbol


# #Test
# # define enhancer region
# enh <- GRanges("chr15:32991002-32994400")
# 
# # extract gene ranges
# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# # txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
# 
# # genes <- genes(txdb, single.strand.genes.only=FALSE)
# genes <- genes(txdb)
#                
# # find overlap
# hits <- findOverlaps(enh, genes)
# gene_ids <- genes[subjectHits(hits)]$gene_id   # Entrez IDs
# print(gene_ids)
# 
# # 4. Map Entrez -> Symbol
# gene_symbols <- mapIds(org.Hs.eg.db,
#                        keys = gene_ids,
#                        column = "SYMBOL",
#                        keytype = "ENTREZID",
#                        multiVals = "first")
# 
# # Show result
# data.frame(EntrezID = gene_ids, Symbol = gene_symbols)


#AllEnhancers
setwd("~/Manuscripts/98RWDisEnhPlus")
Disease2EnhancerFile = "Data/Disease2Enhancers.txt"
Disease2Enhancer = read.delim(Disease2EnhancerFile, sep = "\t", header = FALSE)

dim(Disease2Enhancer)
AllEnhancers = vector()
for(i in 1:nrow(Disease2Enhancer)){
  eidstr = Disease2Enhancer[i,3]
  eidvec = strsplit(eidstr,", ")[[1]]
  for(j in 1:length(eidvec)){
    eid = eidvec[j]
    AllEnhancers = append(AllEnhancers, eid)
  }
}
AllEnhancers = unique(AllEnhancers)

length(AllEnhancers) 

library(hash)

suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))

enh2EntrezID = hash()
for (en in AllEnhancers){
  print(en)
  enh <- GRanges(en)
  
  # extract gene ranges
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  genes <- genes(txdb)
  
  # find overlap
  hits <- findOverlaps(enh, genes)
  gene_ids <- genes[subjectHits(hits)]$gene_id   # Entrez IDs
  if(!is.null(gene_ids)){
    print(gene_ids)
    enh2EntrezID[enh] = gene_ids  
  } 
}
print(length(enh2EntrezID))

# Convert hash to data.frame
keys_vec <- keys(enh2EntrezID)
vals_vec <- values(enh2EntrezID)

df <- data.frame(
  key = keys_vec,
  value = sapply(vals_vec, function(v) paste(v, collapse = ",")),
  stringsAsFactors = FALSE
)

# Write to TXT (tab-separated)
write.table(df, file = "Data/enh2EntrezID.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
