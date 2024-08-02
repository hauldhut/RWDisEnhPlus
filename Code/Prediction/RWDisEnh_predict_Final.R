Method = "RWDisEnhPlus"#RWDisEnh/RWDisEnhPlus

start_time <- Sys.time()

library(RandomWalkRestartMH)
library(igraph)
library(foreach)
library(doParallel)

setwd("~/Manuscripts/98RWDisEnh2/Code/Prediction")

#Shared Gene
enhancer1 <- read.delim("../../Data/EnhNet_SharedGene.txt",header = FALSE)
enhancer1.frame <- data.frame(enhancer1[[1]], enhancer1[[3]])
enhancer1.g <- graph.data.frame(d = enhancer1.frame, directed = FALSE)

#Sequence
enhancer2 <- read.delim("../../Data/EnhNet_Sequence_All.txt",header = FALSE)
enhancer2.frame <- data.frame(enhancer2[[1]], enhancer2[[3]])
enhancer2.g <- graph.data.frame(d = enhancer2.frame, directed = FALSE)
enhancer2.weight = enhancer2[[2]]
E(enhancer2.g)$weight <- enhancer2.weight


if(Method == "RWDisEnhPlus"){
  enhancer_MultiplexObject <- create.multiplex(list(enhancer1.g,enhancer2.g),Layers_Name = c("enhancer1","enhancer2"))  
}else{
  enhancer_MultiplexObject <- create.multiplex(list(enhancer1.g),Layers_Name = c("enhancer"))
}


#Add disease nw
disease <- read.delim("../../Data/DOBasedOMIMEntitySimilarityNet.txt",header = FALSE)

disease.frame <- data.frame(disease[[1]], disease[[3]])
disease.weight = disease[[2]]

disease.g <- graph.data.frame(d = disease.frame, directed = FALSE)
E(disease.g)$weight <- disease.weight

disease_MultiplexObject <- create.multiplex(list(disease.g),Layers_Name = c("disease"))

#Add EDRelation
ED.frame <- read.csv("../../Data/EDRelation.csv", header = TRUE)
ED.frame <- ED.frame[which(ED.frame$enhancer %in% enhancer_MultiplexObject$Pool_of_Nodes),]
ED.frame <- ED.frame[which(ED.frame$disease %in% disease_MultiplexObject$Pool_of_Nodes),]

#Create multiplex-heterosgenous nw
enhancer_Disease_Net <- create.multiplexHet(enhancer_MultiplexObject, disease_MultiplexObject, 
                                            ED.frame)

enhancerHetTranMatrix <- compute.transition.matrix(enhancer_Disease_Net)

#pick unique disease
unique_disease <- unique(ED.frame$disease)

#For a fair comparison with H (on the same set of diseases, only used EnhNet_SharedGene.txt)
enhancer1_MultiplexObject <- create.multiplex(list(enhancer1.g),
                                             Layers_Name = c("enhancer1"))
ED.frame1 <- read.csv("../../Data/EDRelation.csv", header = TRUE)
ED.frame1 <- ED.frame1[which(ED.frame1$enhancer %in% enhancer1_MultiplexObject$Pool_of_Nodes),]
ED.frame1 <- ED.frame1[which(ED.frame1$disease %in% disease_MultiplexObject$Pool_of_Nodes),]
dim(ED.frame1)
unique_disease <- unique(ED.frame1$disease)

#Extract DOID2TraitMap
library(hash)
Disease2EnhancerFile = "../../Data/Disease2Enhancers.txt"
Disease2Enhancer = read.delim(Disease2EnhancerFile, sep = "\t", header = FALSE)
DOID2TraitMap = hash()

for(i in 1:nrow(Disease2Enhancer)){
  did = Disease2Enhancer[i,1]
  trait = Disease2Enhancer[i,2]
  DOID2TraitMap[[did]] = trait
}



#loop through
library(hash)
h.d2ranking = hash()

for(i in 1:length(unique_disease)){ 
  
  library(RandomWalkRestartMH)
  library(igraph)
  library(ROCR)
  
  seeddisease = unique_disease[[i]]
  
  cat("==> Disease ", i, "/", length(unique_disease), ": ", seeddisease,"\n")
  
  disease_relation = ED.frame[which(ED.frame$disease==seeddisease),]
  
  SeedEnhancer = disease_relation$enhancer[c(which(disease_relation$disease==seeddisease))]
  
  #compute 
  RWRH_enhancer_Disease_Results <- Random.Walk.Restart.MultiplexHet(enhancerHetTranMatrix,
                                                                    enhancer_Disease_Net,SeedEnhancer,
                                                                    seeddisease, r = 0.5)
  
  tf = RWRH_enhancer_Disease_Results$RWRMH_Multiplex1
  h.d2ranking[[seeddisease]] = tf  
}

#####################################################
#Summarize for each k
k=10
h.K2TotalEvidencedNewAssoc = hash()
h.K2topKEnhEvidence = hash()

for(k in seq(10, 100, by = 10)){
  # if(k!=10) next
  cat("k =",k,"\n")
  df.topKEnh <- data.frame(matrix(ncol = 2+k, nrow = 0))
  
  for(did in unique_disease){ 
    tf = h.d2ranking[[did]]
    df.topKEnh[nrow(df.topKEnh)+1,] = c(as.character(did),as.character(DOID2TraitMap[[did]]),as.vector(tf$NodeNames[1:k]))
  }
  
  write.table(df.topKEnh, paste0("../../Results/",Method,"_predict_top",k,".txt"), na ="", row.names=FALSE, col.names = FALSE, sep='\t', quote=FALSE)
  
  #Extract Enh2DONameMap, Enh2DOIDMap
  Enh2DONameMap = hash()
  Enh2DOIDMap = hash()
  for(i in 1:nrow(df.topKEnh)){
    did = df.topKEnh[i,1]
    trait = df.topKEnh[i,2]
    for(j in 1:k){
      eid = df.topKEnh[i,2+j]
      traitset = vector()
      didset = vector()
      if(eid %in% keys(Enh2DONameMap)){
        traitset = Enh2DONameMap[[eid]]
        didset = Enh2DOIDMap[[eid]]
      }
      traitset = append(traitset,trait)
      didset = append(didset,did)
      
      traitset = unique(traitset)
      didset = unique(didset)
      
      Enh2DONameMap[[eid]] = traitset
      Enh2DOIDMap[[eid]] = didset
    }
  }
  
  #Extract GWASTrait2DOIDMap
  GWASTrait2DOIDMapFile = "../../Data/TraitSet_2_DOID.txt"
  GWASTrait2DOID = read.delim(GWASTrait2DOIDMapFile, sep = "\t", header = FALSE)
  
  GWASTrait2DOIDMap = hash()
  
  for(i in 1:nrow(GWASTrait2DOID)){
    gwastrait = GWASTrait2DOID[i,1]
    doidvec = strsplit(GWASTrait2DOID[i,3],", ")
    GWASTrait2DOIDMap[[gwastrait]] = doidvec
  }
  
  #findEvidenceForEnhancerFromGWAS_Ver2
  WorkingDir <- "../../Data/AllEnhancers"
  
  df.topKEnhEvidence = NULL
  
  
  TotalHit=0
  TotalEvidencedNewAssoc=0
  EvidenceDOIDSet = vector()
  GWASTraitSet = vector()
  #Traverse topKEnh (Enh2DOIDMap)
  for(eid in keys(Enh2DOIDMap)){
    DOIDSet = Enh2DOIDMap[[eid]]
    EnhDisAssocFile = paste0(WorkingDir,"/",eid,".txt")
    if(!file.exists(EnhDisAssocFile)) next
    
    EnhInfo = read.delim(EnhDisAssocFile, sep = "\t", header = TRUE)
    c = 0
    evidenceddoidset = vector()
    
    # cat(eid, nrow(EnhInfo),"\n")
    
    for(i in 1:nrow(EnhInfo)){
      trait = tolower(EnhInfo$trait[i])
      # cat(trait,"\n")
      if(!(trait %in% keys(GWASTrait2DOIDMap))) next
      GWASTraitSet = append(GWASTraitSet,trait)
      
      GWASDOIDSet = GWASTrait2DOIDMap[[trait]]
      DOID_intersect = intersect(DOIDSet, GWASDOIDSet)
      
      if(length(DOID_intersect)>0){
        # cat(paste(DOID_intersect,collapse = " "),"\n")
        cat(eid,"\t",trait,"\t",paste(DOID_intersect,collapse = " "),"\t",EnhInfo$rsid[i],"(",EnhInfo$p[i],")","\t",EnhInfo$pmid[i],"\t",EnhInfo$consequence[i],"\t",EnhInfo$hgnc[i],"\n")
        
        df.topKEnhEvidence = rbind(df.topKEnhEvidence, data.frame(trait=trait,eid=eid,hgnc=EnhInfo$hgnc[i],rsid = paste0(EnhInfo$rsid[i]," (",EnhInfo$p[i],")"),consequence=EnhInfo$consequence[i],pmid=EnhInfo$pmid[i]))
        
        evidenceddoidset = append(evidenceddoidset,DOID_intersect)
      }
      c=c+length(DOID_intersect)
    }
    evidenceddoidset = unique(evidenceddoidset)
    if (c>0){
      TotalHit=TotalHit+c;
      TotalEvidencedNewAssoc=TotalEvidencedNewAssoc+length(evidenceddoidset);
      EvidenceDOIDSet = append(EvidenceDOIDSet,evidenceddoidset)
    }
  }
  
  
  
  GWASTraitSet = unique(GWASTraitSet)
  EvidenceDOIDSet = unique(EvidenceDOIDSet)
  
  cat("DOID2TraitMap:",length(DOID2TraitMap),"\n")
  cat("Enh2DOIDMap:",length(Enh2DOIDMap),"\n")
  cat(length(GWASTraitSet), "\t", paste(GWASTraitSet,collape=" "),"\n")
  cat("TotalHit:",TotalHit,"\n")
  cat("TotalEvidencedNewAssoc:", TotalEvidencedNewAssoc,"\n")
  cat("EvidenceDOIDSet: ", length(EvidenceDOIDSet),"\t",paste(EvidenceDOIDSet,collape=" "),"\n")
  
  write.csv(df.topKEnhEvidence, paste0("../../Results/",Method,"_predict_top",k,"_Evidence.csv"), na ="", row.names=FALSE, quote=FALSE)
  
  h.K2TotalEvidencedNewAssoc[[as.character(k)]] = TotalEvidencedNewAssoc
  h.K2topKEnhEvidence[[as.character(k)]] = df.topKEnhEvidence
}
h.K2TotalEvidencedNewAssoc

#Summarize Evidence
k=10
for(k in seq(10, 100, by = 10)){
  df.k = h.K2topKEnhEvidence[[as.character(k)]]
  
  h.trait_eid_pmid_hgnc2rsid = hash()
  v.trait_eid_hgnc = vector()
  for(i in 1:nrow(df.k)){
    trait_eid_pmid_hgnc = paste0(df.k$trait[i],"_",df.k$eid[i],"_",df.k$pmid[i],"_",df.k$hgnc[i])
    trait_eid_hgnc = paste0(df.k$trait[i],"_",df.k$eid[i],"_",df.k$hgnc[i])
    v.trait_eid_hgnc = append(v.trait_eid_hgnc, trait_eid_hgnc)
    
    rsid = df.k$rsid[i]
    rsidset = vector()
    if(trait_eid_pmid_hgnc %in% keys(h.trait_eid_pmid_hgnc2rsid)){
      rsidset = h.trait_eid_pmid_hgnc2rsid[[trait_eid_pmid_hgnc]]
    }
    rsidset = append(rsidset,rsid)
    rsidset = unique(rsidset)
    h.trait_eid_pmid_hgnc2rsid[[trait_eid_pmid_hgnc]] =  rsidset
  }
  v.trait_eid_hgnc = unique(v.trait_eid_hgnc)
  cat(k,": NuHit:",length(v.trait_eid_hgnc),"\n")
  
  df.topKEnhEvidenceSum = NULL
  for(trait_eid_pmid_hgnc in keys(h.trait_eid_pmid_hgnc2rsid)){
    s = strsplit(trait_eid_pmid_hgnc,"_")[[1]]
    rsidset = h.trait_eid_pmid_hgnc2rsid[[trait_eid_pmid_hgnc]]
    rsidlist = paste(rsidset, collapse = ", ")
    df.topKEnhEvidenceSum = rbind(df.topKEnhEvidenceSum, data.frame(trait = s[1], eid = s[2],hgnc=s[4], rsid=rsidlist, pmid=s[3]))
  }
  df.topKEnhEvidenceSum
  
  write.table(df.topKEnhEvidenceSum, paste0("../../Results/",Method,"_predict_top",k,"_EvidenceSum.txt"), na ="", row.names=FALSE, quote=FALSE, sep="\t")
}
