Method = "RWDisEnh"#RWDisEnh/RWDisEnhPlus

start_time <- Sys.time()

library(RandomWalkRestartMH)
library(igraph)

#need to install foreach and doParallel packages for this code to run
library(foreach)
library(doParallel)

setwd("~/Manuscripts/98RWDisEnhPlus/Code/K-Fold")

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

#add func for RWR on multiplex-heter nw
do_something <- function(enhancer_MultiplexObject,disease_MultiplexObject,
                         EDRelation,SeedEnhancer, seeddisease, prd_enhancers) {
  
  #Create multiplex-heterosgenous nw
  
  EDRelation_enhancer <- EDRelation[which(EDRelation$enhancer %in% enhancer_MultiplexObject$Pool_of_Nodes),]
  
  #Create multiplex-heterosgenous nw
  enhancer_Disease_Net <- create.multiplexHet(enhancer_MultiplexObject, disease_MultiplexObject, 
                                              EDRelation_enhancer)
  
  enhancerHetTranMatrix <- compute.transition.matrix(enhancer_Disease_Net)
  
  #compute 
  # tau <- c(1,1)
  RWRH_enhancer_Disease_Results <- Random.Walk.Restart.MultiplexHet(enhancerHetTranMatrix,
                                                                    enhancer_Disease_Net,SeedEnhancer,
                                                                    seeddisease, r = 0.5)
  
  #create labels for enhancer results
  # tf = RWRH_enhancer_Disease_Results$RWRMH_Results_MultiplexNodes
  tf = RWRH_enhancer_Disease_Results$RWRMH_Multiplex1
    
  tf$labels <- ifelse(tf$NodeNames %in% prd_enhancers, 1, 0)
  
  # calculating AUC for each enhancer of a disease
  tf$Score[is.na(tf$Score)] = 0
  
  resultspred = prediction(tf$Score, tf$labels)
  
  pauc.perf = performance(resultspred, measure = "auc")
  
  return(pauc.perf@y.values[[1]])
}

#count enhancer for each disease
sub_sum <- aggregate(enhancer~disease, data=ED.frame, FUN=function(x) c(count=length(x)))

#extract disease with only 3 or more enhancer(k=3)
k=3
sub_sum <- sub_sum[which(sub_sum$enhancer>=k),]
sub_sum$disease_no <- c(1:length(sub_sum$disease))

#extract ED.frame with only disease from sub_sum
ED.frame1 <- ED.frame[which(ED.frame$disease %in% sub_sum$disease),]
rownames(ED.frame1) <- NULL #reset frame index

#func to assign k groups for each set of disease-enhancer, 
#as well as increment group no. for each group (k=3)
assign_group_no <- function(sub_sum,ED.frame1,k) {
  
  #set an empty data frame for a new ED.frame
  mylist.names <- c("disease","enhancer", "disease_no","group_no")
  ED.frame2 <- sapply(mylist.names,function(x) NULL)
  
  
  for (j in 1:length(sub_sum$disease)) {
    
    count = sub_sum$enhancer[[j]]
    
    if(count<k) next
    
    set_no = floor(count/k)
    
    print(paste(k,count, set_no))
    
    group_vec = vector()
    for(gi  in 1:(k-1)){
      group_vec = c(group_vec,rep(gi,set_no))
    }
    group_vec = c(group_vec,rep(k,count-set_no*(k-1)))

    # group_vec <- c(rep(1,set_no), rep(2,set_no), rep(3,set_no), rep(4,set_no),rep(5,set_no), rep(6,set_no),rep(7,set_no), rep(8,set_no), rep(9,set_no),rep(10,(count-set_no*9)))
    
    subset <- ED.frame[which(ED.frame$disease==sub_sum$disease[[j]]),]
    subset$disease_no <- rep(j,count)
    subset$group <- group_vec
    
    ED.frame2 <- rbind(ED.frame2,subset)
  }
  return(ED.frame2)
}

#assign each group of n/k enhancer-disease with an group id
out <- assign_group_no(sub_sum,ED.frame1,k)
out$group <- (out$disease_no-1)*k+out$group
ED.frame2 <- out

#set an empty data frame for a new ED.frame
auc_results <- sapply(c("disease","group","auc"),function(x) NULL)

for (i in 1:max(ED.frame2$group)) {
  seeddisease = unique(ED.frame2$disease[which(ED.frame2$group==i)])
  auc_results$disease[i] <- seeddisease
  auc_results$group[i] <- i
}

#set up paralell processing (adjust the no_cores as per running system)
no_cores <- 6
cl <- makeCluster(no_cores)
registerDoParallel(cl)

#loop through
auc <- foreach(i = 1:max(ED.frame2$group), .combine = rbind) %dopar% {
  
  library(RandomWalkRestartMH)
  library(igraph)
  library(ROCR)
  
  prd_enhancers = ED.frame2$enhancer[which(ED.frame2$group==i)]
  seeddisease = unique(ED.frame2$disease[which(ED.frame2$group==i)])
  
  disease_relation = ED.frame2[which(ED.frame2$disease==seeddisease),]
  SeedEnhancer = disease_relation$enhancer[-c(which(disease_relation$enhancer %in% prd_enhancers))]
  
  # get bipartite graph without prd_enhancers - disease linkages
  EDRelation <- ED.frame2[-with(ED.frame2, which(enhancer %in% prd_enhancers & disease %in% seeddisease)),][1:2]
  
  n <- do_something(enhancer_MultiplexObject,disease_MultiplexObject,
                    EDRelation,SeedEnhancer, seeddisease, prd_enhancers)
}

stopCluster(cl)

a <- split(auc, rep(1:ncol(auc), each = nrow(auc)))
auc_results$auc <- a$`1`
write.csv(auc_results, paste0("../../Results/",Method,"_kfold.csv"))

sum(auc_results$auc)/length(auc_results$auc)

#group and calculate by disease
final <- aggregate(auc~disease, data=auc_results, FUN=function(x) c(mean=mean(x), count=length(x)))
write.csv(final, paste0("../../Results/",Method,"_kfold_byDisease.csv"))



end_time <- Sys.time()
end_time - start_time