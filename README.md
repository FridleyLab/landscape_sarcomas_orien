# Genomic, transcriptomic, and immunogenomic landscape of over 1,300 sarcomas of diverse histology subtypes: A report from the ORIEN collaborative network

This repository contains code to reproduce the results of the manuscript "Genomic, transcriptomic, and immunogenomic landscape of over 1,300 sarcomas of diverse histology subtypes: a report from the ORIEN collaborative network". The manuscript can be obtained from [LINK/DOI].

Data used in this code is not contained in this repository and is managed by the ORIEN network. When feasible, we have added details about the data formats to help the reader understand the code. Analysis of copy number variation and mutational data was carried out by [Dr. Alex Soupir](https://www.alexsoupir.com/). Analyses of gene expression data was carried out by Oscar Ospina.

The .fastq files from Illumina outputs were processed using Space Ranger. The Seurat package was used to import and further process gene expression counts. Biological identification of the spots was achieved using STdeconvolve. Gene expression gradients (STgradient) and spatial gene set enrichment analysis (STenrich) were conducted in spatialGE.

## `gene_expression_analyses`
The `Gene expression analyses` folder contains code to perform exploratory data analysis, unsupervised clustering, differential gene expression, and gene set enrichment analysis of the sarcoma RNAseq samples in the cohort. Code is also provided for immune cell type deconvolution and survival analysis using Cox models.
* `gene_expression_analyses/sarcoma_avatar_rnaseq_data_exploration_gene_level.Rmd`: Code to assess number and distribution of samples in the cohort.
* `gene_expression_analyses/sarcoma_avatar_rnaseq_mcpdeconv`: Code to estimate relative abundance of immune cell types using MCPdeconv RNAseq deconvolution.
