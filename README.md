# RWDisEnh+
Enhancers regulate the transcription of target genes, influencing their expression levels. Consequently, aberrant enhancer activity can lead to diseases. To date, millions of enhancers have been identified, yet only a small fraction of them have been linked to diseases. This underscores the pressing need to develop computational methods to predict associations between diseases and enhancers. In a previous study, we hypothesized that enhancers sharing target genes could be associated with similar diseases; we thus proposed a network-based method utilizing a heterogeneous network that connected a shared gene-based enhancer network and a disease similarity network to predict novel disease-associated enhancers. 

In this study, we extended our hypothesis, proposing that similar enhancers could be associated with similar diseases to improve prediction accuracy. We additionally built a sequence-based enhancer similarity network and integrated it with the existing heterogeneous network to form a multiplex-heterogeneous network of diseases and enhancers. Furthermore, our random walk with restart scheme on the heterogeneous network (RWDisEnh) was extended to operate on the multiplex-heterogeneous network (RWDisEnh+) to measure the degree of association. Experimental results show that RWDisEnh+ outperforms RWDisEnh in terms of AUC values obtained by 3-fold cross-validation; as well as the number of top-ranked enhancers supported by direct evidence. Our findings demonstrate that RWDisEnh+ can predict novel enhancers associated with specific diseases.

## Construction of networks of diseases and enhancers
![Construction of networks of diseases and enhancers](https://github.com/hauldhut/RWDisEnh/blob/main/Figure1.png)

## Data
* **Disease2Enhancers.txt**:
  - All disease-enhancer associations collected from DiseaseEnhancer database, in which disease names were mapped to DO (Disease Ontology) ID
* **EDRelation.csv**:
  - All binary disease-enhancer association collected from DiseaseEnhancer database
* **DOBasedOMIMEntitySimilarityNet.txt**:
  - DO-based disease similarity network.
* **EnhNet_SharedGene.txt**:
  - Shared gene-based enhancer network
* **EnhNet_Sequence_All.txt**:
  - Sequence-based enhancer similarity network
* **Trait_2_DOID.txt**:
  - A mapping between disease (trait) to Disease Ontology identifier (DOID), used to collect evidence for top ranked enhancers for each disease
* **AllEnhancers.zip**:
  - The file contains enhancer-phenotype associations collected from GWASCatalog using PhenoScanner for all enhancers in DiseaseEnhancer database 

## Code
* Prepare
  - **Install R packages**: *RandomWalkRestartMH, igraph, foreach, doParallel, ROCR, Metrics, hash, phenoscanner, clusterProfiler, org.Hs.eg.db*

## Results
* Contains all outputs/results of k-fold cross-validation and prediction


