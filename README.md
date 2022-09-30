# RWDisEnh
Motivation: Enhancers regulate the transcription of target genes, causing a change in expression level. Thus, the aberrant activity of enhancers can lead to diseases. To date, millions of enhancers have been identified, yet a small portion of them have been found to be associated with diseases. This raises a pressing need to develop computational methods to predict associations between diseases and enhancers. In a previous study, we assumed that enhancers sharing target genes could be associated with similar diseases; thus, proposed a network-based method based on a heterogeneous network connecting a shared gene-based enhancer network and a disease similarity network to predict novel disease-associated enhancers.
Results: In this study, we assumed that similar enhancers could be associated with similar diseases to predict the association. Thus, we additionally built a sequence-based enhancer similarity network and integrated it with the heterogeneous network to form a multiplex-heterogeneous network of diseases and enhancers. In addition, the random walk with restart scheme on the heterogeneous, RWDisEnh, is extended work on the multiplex-heterogeneous network, RWDisEnh+, to measure the degree of the association. Experimental results show that RWDisEnh+ outperforms RWDisEnh in terms of AUC values obtained by 3-fold cross-validation as well as the number of top-ranked enhancers supported by direct evidence.![image](https://user-images.githubusercontent.com/17016237/193174437-19b371b2-8a9c-42bf-94b9-08218ac4f707.png)

# Data
## Disease2Enhancers.txt: All disease-enhancer associations collected from DiseaseEnhancer database, in which disease names were mapped to DO (Disease Ontology) ID
## EDRelation.csv: All binary disease-enhancer association collected from DiseaseEnhancer database
## DOBasedOMIMEntitySimilarityNet.txt: DO-based disease similarity network.
## EnhNet_SharedGene.txt: Shared gene-based enhancer network
## EnhNet_Sequence_All.txt: Sequence-based enhancer similarity network
## AllEnhancers: Contains enhancer-phenotype associations collected from GWASCatalog using PhenoScanner for all enhancers in DiseaseEnhancer database 

# Code
## Data Preprocessing
### QueryFromGWASCatalog_ByEnhancer.R: Collect enhancer-phenotype associations from GWASCatalog using PhenoScanner for all enhancers in DiseaseEnhancer database
## KFold
### 
## Prediction
### RWDisEnh_predict.R: Predicts, selects top k-ranked enhancers, collects and summarizes evidence for the top-k ranked enhancers for the heterogeneous and multiplex-heterogeneous networks



# Results
## Contains all outputs/results of k-fold validation/prediction


