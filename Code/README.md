# Code
* Data Preprocessing
  - **QueryFromGWASCatalog_ByEnhancer.R**: Collect enhancer-phenotype associations from GWASCatalog using PhenoScanner for all enhancers in DiseaseEnhancer database
    
* KFold cross-validation
  - **RWDisEnhPlus_kfold_Final.R**: Performs 3-fold cross-validation for the heterogeneous and multiplex-heterogeneous networks  
    - Change Method to RWDisEnh or RWDisEnhPlus for RWDisEnh or RWDisEnh+, respectively

* Prediction
  - **RWDisEnhPlus_predict_Final.R**: Predicts, selects top k-ranked enhancers, collects and summarizes evidence for the top-k ranked enhancers for the heterogeneous and multiplex-heterogeneous networks
    - Change Method to RWDisEnh or RWDisEnhPlus for RWDisEnh or RWDisEnh+, respectively
      
* Analysis
  - **RWDisEnhPlus_Summarize_topkEvidence.R**: Summarize the number of evidence found for top-ranked enhancers.
  - **RWDisEnhPlus_Enrich_Pathways_Final.R**: Perform KEGG Pathway Enrichment Analysis for involved SNPs

