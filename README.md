# Ayesha-Abdelazim_5_Assignment
E-GEOD-36980: Alzheimer's Disease Brain Differential Gene Expression Analysis
Microarray analysis of post-mortem Alzheimer's disease brain samples (79 total: 47 Control vs 32 Alzheimer's) using Bioconductor in R.

## Analysis Workflow
- **Quality Control**: arrayQualityMetrics for pre/post-normalization QC
- **Normalization**: RMA (Robust Multi-array Average)
- **Filtering**: Removal of probes with low intensity across samples
- **Differential Expression**: Limma with empirical Bayes moderation
- **Visualization**: ggplot2 for volcano plots, pheatmap for heatmaps
