# BIOS841_Statistical_Consulting_Final_Project

This repository contains the data, analyses, and final report for a statistical consulting class.

## Detecting changes in splicing efficiency across the mouse genome
Principle Investigator: Dr. Mauro Calabrese
Coauthors: Renee Ge, Xueyao Wang, Laura Watson, Cuiran Shi, Longfei Zhang

This study investigates patterns of RNA splicing efficiency across various tissues in the mouse genome. The first aim is to identify genes that exhibit significant changes in splicing efficiency across a collection of different tissues. The second aim is to explore whether there are discernible patterns or correlations among the genes shown to have these significant changes, potentially revealing broader regulatory or functional trends in splicing across tissues.

## Data 

The RNA sequencing data utilized in this study were derived from the dataset associated with the publication “Mapping the mouse Allelome reveals tissue-specific regulation of allelic expression,” submitted to NCBI (GEO). This dataset was further enriched with additional RNA-seq data provided by the Calabrese Lab at UNC Chapel Hill

## Methods

For our first analysis, we performed differential splicing analysis for all genes across mouse tissues at different developmental stages. Our second analysis addresses the question of determining if there are patterns or correlations among the genes with significant changes in splicing efficiency. We applied Weighted Gene Co-expression Network Analysis (WGCNA) to identify modules of genes with coordinated splicing efficiency patterns across samples.

## Project Organization
* 25_03_16_Parsed_Revised_full.csv.zip is a compressed file containing the raw data file
* data_processing.R takes in 25_03_16_Parsed_Revised_full.csv and outputs a processed data file "processed_data.csv" to be used in analysis
* processed_data_NAs_included.csv is the same as the processed_data.csv file but the missing values are not excluded
* DESeq2_Analysis.R contains R code that conducts differential splicing analysis for all tissues included in the processed data file and saves all output to the DESeq_Results folder
