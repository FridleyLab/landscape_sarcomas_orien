# Genomic, transcriptomic, and immunogenomic landscape of over 1,300 sarcomas of diverse histology subtypes: A report from the ORIEN collaborative network

This repository contains code to reproduce the results of the manuscript "Genomic, transcriptomic, and immunogenomic landscape of over 1,300 sarcomas of diverse histology subtypes: a report from the ORIEN collaborative network". The manuscript can be obtained from [LINK/DOI].

Data used in this code is not contained in this repository and is managed by the ORIEN network. When feasible, we have added details about the data formats to help the reader understand the code. Analysis of copy number variation and mutational data was carried out by [Dr. Alex Soupir](https://www.alexsoupir.com/). Analyses of gene expression data were carried out by [Dr. Oscar Ospina](https://github.com/oospina).

## [`somatic_mutations`](https://github.com/FridleyLab/landscape_sarcomas_orien/tree/main/somatic_mutations)
Processing of somatic mutation data from the ORIEN network's sarcoma cohort can be found in the `somatic_mutations` folder. This includes the filtering of VCF somatic mutation files, annotation with AnnoVar, and then comparision/plots of the resulting somatic mutations.

## [`copy_number`](https://github.com/FridleyLab/landscape_sarcomas_orien/tree/main/copy_number)
Scripts for using both gene- and segment-level copy number can be found in the `copy_number` folder. Code at the bottom of `5.0.Comparing_total_perMB.R` is modified from the `compareTCGA()` function within the *maftools* package to handle samples with differing capture technologies/lengths.

## [`fusions`](https://github.com/FridleyLab/landscape_sarcomas_orien/tree/main/fusions)
Fusion data was processed after all other scripts were written then used to update the clinical information with histology reassignments. After updating, all scripts were rerun with the update histology subtypes. 

## [`gene_expression_analyses`](https://github.com/FridleyLab/landscape_sarcomas_orien/tree/main/gene_expression_analyses_code)
The `Gene expression analyses` folder contains code to perform exploratory data analysis, unsupervised clustering, differential gene expression, and gene set enrichment analysis of the sarcoma RNAseq samples in the cohort. Code is also provided for immune cell type deconvolution and survival analysis using Cox models.
* `gene_expression_analyses_code/sarcoma_manuscript_figures.Rmd`: Code to to reproduce the figures pertaining to RNAseq analysis. The script takes outputs from other scripts included here.
