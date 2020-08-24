# Installation

### Environment & Dependencies

- [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest) V. 3.0.2

This pipeline analysis use R studio V. 4.0.2 and the packages: 

- [SoupX](https://github.com/constantAmateur/SoupX) V. 1.4.5
- [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder) V. 2.0.3
- [Seurat](https://satijalab.org/seurat/) V. 3.1.2
- [Monocle3](https://cole-trapnell-lab.github.io/monocle3/docs/installation/) V. 0.2.3.0
- [gProfiler2](https://biit.cs.ut.ee/gprofiler/page/r) V. 0.1.9



# Data

Raw files are deposited in ArrayExpress, under accession number E-MTAB-9467



# Data Mapping 

[Docker](https://hub.docker.com/r/jantarika/cellranger_denguetimecourse) is provided

[Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest) V. 3.0.2 and human reference GRCh38-1.2.0 are applied 

Running command as follows: 

cellranger count --id= sample id \
                              --transcriptome=PATH_TO_REF \
                              --fastqs=PATH_TO_RAW_FILES \
                              --sample=sample name \
                              --chemistry=SC3Pv2 \
                              --expect-cells=5000



# Workflow of downstream analyses

[Docker](https://hub.docker.com/r/jantarika/rstudio_denguetimecourse) is provided

1. Quality control 
   - Correct the expression profiles using SoupX
   - Exclude the cells that have high percentage of mitochondrial counts
   - Exclude doublet using DoubletFinder  

2. Data integration

   - Integrate all the samples using SCTransform normalization generated by Seurat

3. Highly variable genes and pathway analyses

   - Construct Pearson Correlation Coefficient and PCAs
   - Construct unsupervised hierarchical clustering from top 500 genes on PCs
   - Investigate GO enrichments of each cluster using gProfiler2 

4. Trajectory and pseudotime analyses 

   - Create Monocle's object using the information from Seurat's object 
   - Construct trajectory and pseudotime using Monocle3
   - Set Naive at the root of pseudotime for CD4+ T cells and B cell subpopulations
   - Set Effector CD8-3 at the root for effector CD8+ subpopulations





















