#!/usr/bin/env Rscript
# Download libraries
library(Seurat)
library(DoubletFinder)
library(stringr)
library(dplyr)
library(tidyverse)


### READ ME ####
# This script is for excluding the low-quality cells
# Corrected expression matrix from SoupX is used as the input
# To get corrected expression matrix, user can be either run "SoupX.R" script or download the files in A-XXXX-n

# Function for running the standard pre-processing on Seurat object
Standard_PreProcessing_QC <- function(input){
  # Calculate percentage of mitochondrial counts
  print("Calculate % of mitochondrial counts")
  input <- PercentageFeatureSet(input, pattern = "^MT-" ,  col.name = "percent.mt")
  
  # Apply sctransform normalization
  print("SCTransform Normalization")
  input <- SCTransform(input, verbose = FALSE)
  
  # Perform linear dimensional reduction
  print("Perform linear dimensional reduction")
  input <- RunPCA(input ,  verbose = FALSE)
  
  # Cluster the cells
  print("Cluster the cells")
  input <- FindNeighbors(input, dims = 1:30)
  input <- FindClusters(input,  verbose = FALSE)
  
  # Run non-linear dimensional reduction (UMAP)
  print("Run non-linear dimensional reduction (UMAP)")
  input <- RunUMAP(input, dims = 1:30 ,  verbose = FALSE)
  return(input)
}

# Create Seurat object
sc_Seurat <- Read10X(data.dir = "/PATH_TO_DIRECTORY/") # Directory of filtered_feature_bc_matrix by CellRanger
sc_Seurat <- CreateSeuratObject(counts = sc_Seurat)
sc_Seurat_procesing <- sc_Seurat

# Run standard pre-processing on Seurat object
sc_Seurat_procesing <- Standard_PreProcessing_QC(sc_Seurat_procesing)

# Exclude cells that have more than 10% of mitochondrial counts
sc_Seurat_procesing <- subset(sc_Seurat_procesing, subset =  percent.mt < 10)

# Extract the cell barcodes that pass % mito cut-off
sc_Seurat_procesing <- subset(sc_Seurat , cells = colnames(sc_Seurat_procesing))

# Run standard pre-processing on Seurat object that contain only the cells that pass % mito cut-off 
sc_Seurat_procesing <- Standard_PreProcessing_QC(sc_Seurat_procesing)

# Run DoubletFinder, predict and exclude doublets 
sweep.res.list <- paramSweep_v3(sc_Seurat_procesing, PCs = 1:30, sct = TRUE,num.cores = 25)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
pK_value <- find.pK(sweep.stats)

# - Extract the pK value at maximum mean-variance-normalized bimodality coefficient (MeanBC)  
pK_value<-pK_value%>%filter(MeanBC==max(MeanBC))%>%dplyr::select(pK)%>%unlist%>%as.character()%>%as.numeric()

# - Predict doublets 
annotations <- sc_Seurat_procesing@meta.data$seurat_clusters
annotations <- as.character(annotations)
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.039*length(sc_Seurat_procesing%>%colnames))  # Assuming 3.9% doublet formation rate according to the 10x protocol
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
sc_Seurat_procesing <- doubletFinder_v3(sc_Seurat_procesing , PCs = 1:30, pN = 0.25, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
colnames(sc_Seurat_procesing@meta.data)[length(sc_Seurat_procesing@meta.data)] <- "DoubletFinder"

# Visualise doublets 
DimPlot(sc_Seurat_procesing , group.by = "DoubletFinder" , cols = c("red" , "black"))

# Exclude doublets
sc_Seurat <- subset(sc_Seurat , cells =  rownames(sc_Seurat_procesing@meta.data[sc_Seurat_procesing@meta.data$DoubletFinder == "Singlet",])) 

# Add information of sample into metadata
sc_Seurat$time <- "time" # Time-course of DENV infection eg "Day-2"
sc_Seurat$severity <- "severity" # Severity of patient eg "DF"
sc_Seurat$severity_time <- paste(sc_Seurat$severity, sc_Seurat$severity)
  
# Save object that already passed QC steps
saveRDS(sc_Seurat , file = "/PATH_TO_DIRECTORY/OBJECT_NAME.rds")








