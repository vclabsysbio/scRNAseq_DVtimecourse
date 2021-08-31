
# scRNAseq_DENVtimecourse

## Description  

This repository contains all required code to reproduce the analyses in our publication.


## Data

Data was deposited in ArrayExpress, under accession number [E-MTAB-9467](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-9467/)


## Raw data processing 

[Cell Ranger (v.3.0.2)](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest) and the reference human genome GRCh38 1.2.0 were applied for read mapping and UMI qualification

```bash
cellranger count --id=sample id \
                 --transcriptome=PATH_TO_REF \
                 --fastqs=PATH_TO_RAW_FILES \
                 --sample=sample name \
                 --chemistry=SC3Pv2 \
                 --expect-cells=5000
```

## R scripts
R sctipts for data analyses in this pubication including
  - Part 01 : Pre-procseeing and quality controls
              We processed data from cellranger filtered quantification matrix and applied the standard pipeline from [Seurat (v.3.1.2)](https://satijalab.org/seurat/) for data normalisation, clustering and dimensionality reduction. To remove the potential contamination of ambient RNAs and doublets, [SoupX (v.1.4.5)](https://github.com/constantAmateur/SoupX) and [DoubletFinder (v.2.0.3)](https://github.com/chris-mcginnis-ucsf/DoubletFinder) were performed to each sample before data integration.
  - Part 02 : Data normalisation and integration
              The SCTransform from [Seurat (v.3.1.2)](https://satijalab.org/seurat/) was used for normalisation prior to integration. 
  - Part 03 : Highly variable genes (“HVGs”) and Biological Process (BP) analyses
  - Part 04 : Trajectory and pseudotime analyses
    


## Downstream analyses

1. Quality control 
   - Correct the expression profiles using SoupX
   - Exclude the cells that have high percentage of mitochondrial counts
   - Exclude doublet using DoubletFinder  

2. Data integration

   - Integrate all the samples using SCTransform normalisation generated by Seurat

3. Highly variable genes and pathway analyses

   - Construct Pearson Correlation Coefficient and PCAs
   - Construct unsupervised hierarchical clustering from top 500 genes on PCs
   - Investigate GO enrichments of each cluster using gProfiler2 

4. Trajectory and pseudotime analyses 

   - Create Monocle's object using the information from Seurat's object 
   - Set Naive at the root of pseudotime for CD4+ T and B cell subpopulations
   - Set Effector CD8-3 at the root for effector CD8+ T cell subpopulations


## Installation

## Environments and Dependencies

- [Cell Ranger (v.3.0.2)](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest) V. 3.0.2 ([Docker](https://hub.docker.com/r/jantarika/cellranger_denguetimecourse) is available) 

R studio v.4.0.2 and the packages ([Docker](https://hub.docker.com/r/jantarika/rstudio_denguetimecourse) is available): 

- [SoupX](https://github.com/constantAmateur/SoupX) V. 1.4.5
- [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder) V. 2.0.3
- [Seurat](https://satijalab.org/seurat/) V. 3.1.2
- [Monocle3](https://cole-trapnell-lab.github.io/monocle3/docs/installation/) V. 0.2.3.0
- [gProfiler2](https://biit.cs.ut.ee/gprofiler/page/r) V. 0.1.9

Data analyses were performed on Ubuntu 16.04.6 LTS






























