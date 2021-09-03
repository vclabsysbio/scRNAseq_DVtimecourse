
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
R scripts for data analysis in this publication including
  - **Part01** : Pre-processing steps and quality control of each sample  
                 The standard pipeline from [Seurat (v.3.1.2)](https://satijalab.org/seurat/) was applied for data normalisation, clustering and dimensionality reduction. To remove the potential contamination of ambient RNAs and doublets, [SoupX (v.1.4.5)](https://github.com/constantAmateur/SoupX) and [DoubletFinder (v.2.0.3)](https://github.com/chris-mcginnis-ucsf/DoubletFinder) were performed before integrating the data.
                 
  - **Part02** : Data normalisation and integration (Figure 1B)  
                 The SCTransform from [Seurat (v.3.1.2)](https://satijalab.org/seurat/) was used for normalisation prior to integration. 
                 
  - **Part03** : Highly variable genes (HVGs) and biological process (BP) analyses (Figure 2A-C)  
                 The Principal Component Analysis (PCA) was performed (Figure 2A). The union of the top 500 genes from PC1 and PC2 (HVGs) were then used to construct the heatmap (Figure 2B). The HVG and BP analyses of each immune cell type (Figure 2C) were used the same Rscript, except PCA and HVGs were constructed using the integrated object from each cell type.
                 
  - **Part04** : Trajectory and pseudotime analyses (Figure 3E)
    

## Environments and Dependencies

- [Cell Ranger (v.3.0.2)](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest); [Docker](https://hub.docker.com/r/jantarika/cellranger_denguetimecourse) is available


R studio v.4.0.2 and the packages used in this publication; [Docker](https://hub.docker.com/r/jantarika/rstudio_denguetimecourse) is available

- [SoupX (v.1.4.5)](https://github.com/constantAmateur/SoupX) 
- [DoubletFinder (v.2.0.3)](https://github.com/chris-mcginnis-ucsf/DoubletFinder) 
- [Seurat (v.3.1.2)](https://satijalab.org/seurat/) 
- [Monocle3 (v.0.2.3.0)](https://cole-trapnell-lab.github.io/monocle3/docs/installation/) 
- [gProfiler2 (v.0.1.9)](https://biit.cs.ut.ee/gprofiler/page/r) 

Data analyses were performed on Ubuntu 16.04.6 LTS


