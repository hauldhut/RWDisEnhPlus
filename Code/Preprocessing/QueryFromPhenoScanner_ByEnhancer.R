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

AllEnhancersFolder <- "../../Data/AllEnhancers"

library('phenoscanner')

setwd(AllEnhancersFolder)
for(i in 1:length(AllEnhancers)){
  eid <-AllEnhancers[i]
  cat(i,"/",length(AllEnhancers),":",eid,"\n")
  res <- phenoscanner(regionquery=eid)
  
  if(nrow(res$results)==0) next()
  fn <-gsub(":","&",eid)  
  write.table(res$results,paste0(eid,".txt"),sep = "\t", quote = FALSE)
}
