# RWDisEnh+
Enhancers are critical regulatory elements that, when disrupted, contribute to disease pathogenesis by altering gene expression. Despite the identification of millions of enhancers, their disease associations remain largely unexplored, necessitating advanced computational methods. In our previous work, we developed RWDisEnh, a network-driven approach that integrates a shared gene-derived enhancer network with a disease similarity network within a heterogeneous framework to predict disease-enhancer associations. 
In this study, we introduce RWDisEnh+, an improved method that enhances prediction accuracy by incorporating a sequence-based enhancer similarity network into a multiplex-heterogeneous framework. Using an extended random walk with restart (RWR) algorithm, RWDisEnh+ leverages the multiplex-heterogeneous network to rank candidate enhancers based on their association with diseases of interest. Experimental results from 3-fold cross-validation demonstrate that RWDisEnh+ achieves an AUC of 0.874, outperforming RWDisEnh’s AUC of 0.819. Furthermore, RWDisEnh+ identifies more evidence-supported enhancers across top-k rankings (k = 10 to 100), successfully predicting novel associations between 10 enhancers and seven diseases, including asthma, rheumatoid arthritis, and type 2 diabetes. Pathway enrichment analysis reveals that these associations are linked to immune, inflammatory, and metabolic pathways, underscoring their biological relevance. These findings highlight RWDisEnh+’s potential to uncover novel disease-enhancer relationships, advancing our understanding of complex disease mechanisms.

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



