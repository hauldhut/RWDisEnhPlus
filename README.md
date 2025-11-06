# RWDisEnh+
Enhancers are critical regulatory DNA elements that, when dysregulated, can contribute to disease pathogenesis by altering gene expression. Although millions of enhancers have been identified through large-scale genomic projects, their associations with human diseases remain largely uncharacterized, emphasizing the need for robust computational approaches. In our previous work, we developed RWDisEnh, a network-based method that integrates a shared gene–based enhancer network (gEnhNet) with a disease similarity network within a heterogeneous framework to predict disease–enhancer associations. 

In this study, we present RWDisEnh+, an enhanced version of RWDisEnh that incorporates a sequence-based enhancer similarity network (sEnhNet) into a multiplex-heterogeneous network to improve prediction performance. Using an extended random walk with restart (RWR) algorithm, RWDisEnh+ allows information to propagate across disease and enhancer layers, leveraging both gene-based and sequence-based similarity features to rank candidate enhancers for each disease. 

Comprehensive evaluations using 3-fold cross-validation demonstrate that RWDisEnh+ achieves an average AUC of 0.874, outperforming RWDisEnh’s AUC of 0.819. Moreover, RWDisEnh+ identifies a larger number of evidence-supported disease–enhancer associations across top-k rankings, including 10 enhancers linked to seven diseases such as asthma, rheumatoid arthritis, and type 2 diabetes. GWAS validation and pathway enrichment analyses further reveal that these predicted associations are enriched in immune, inflammatory, and metabolic pathways, highlighting their biological relevance. Overall, RWDisEnh+ provides a stable and effective framework for predicting novel disease–enhancer associations, offering new insights into enhancer-mediated gene regulation and the genetic architecture of complex diseases.

## Construction of networks of diseases and enhancers
![Construction of networks of diseases and enhancers](https://github.com/hauldhut/RWDisEnh/blob/main/Figure1.png)

## Repo structure
- **Data**: Contains all data 
- **Code**: Contains all source code to reproduce all the results
- **Results**: Contains all outputs/results of k-fold cross-validation and prediction

## How to run
- Install R packages
  - *RandomWalkRestartMH, igraph, foreach, doParallel, ROCR, Metrics, hash, phenoscanner, biomaRt, clusterProfiler, org.Hs.eg.db*
- Download the repo
- Follow instructions in the folder **Code** to run



