
# scRNAseq_DENVtimecourse

## Description  

This repository contains the required codes to perform scRNA-sequencing analysis in our publication entitled "Single-cell Temporal Analysis of Natural Human Dengue Virus Infection Reveals Expansion of Skin Homing Lymphocyte Subsets One Day before Defervescence"


## Data

Data was deposited in the ArrayExpress repository under the accession number [E-MTAB-9467](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-9467/)


## Raw data processing 

Sequenced data was analysed using [CellRanger (v.3.0.2)](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest) and reference human genome GRCh38 1.2.0

```bash
cellranger count --id = <sample_id> \
                 --transcriptome = <refdata-cellranger-GRCh38-1.2.0> \
                 --fastqs = <fastq_path> \
                 --sample = <sample_name> \
                 --chemistry = SC3Pv2 \
                 --expect-cells = 5000
```

## R scripts
R sctipts for data analysis in this pubication including
  - Part 01 : Pre-procseeing and quality controls
  We processed data from cellranger filtered quantification matrix and applied the standard pipeline from [Seurat (v.3.1.2)](https://satijalab.org/seurat/) for data normalisation, clustering and dimensionality reduction. To remove the potential contamination of ambient RNAs and doublets, [SoupX (v.1.4.5)](https://github.com/constantAmateur/SoupX) and [DoubletFinder (v.2.0.3)](https://github.com/chris-mcginnis-ucsf/DoubletFinder) were performed in each sample before data integration.
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



## Environments and Dependencies

- [Cell Ranger (v.3.0.2)](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest) : [Docker](https://hub.docker.com/r/jantarika/cellranger_denguetimecourse) is available


R studio v.4.0.2 and the packages : [Docker](https://hub.docker.com/r/jantarika/rstudio_denguetimecourse) is available

- [SoupX (v.1.4.5)](https://github.com/constantAmateur/SoupX) 
- [DoubletFinder (v.2.0.3)](https://github.com/chris-mcginnis-ucsf/DoubletFinder) 
- [Seurat (v.3.1.2)](https://satijalab.org/seurat/) 
- [Monocle3 (v.0.2.3.0)](https://cole-trapnell-lab.github.io/monocle3/docs/installation/) 
- [gProfiler2 (v.0.1.9)](https://biit.cs.ut.ee/gprofiler/page/r) 

Data analyses were performed on Ubuntu 16.04.6 LTS






























