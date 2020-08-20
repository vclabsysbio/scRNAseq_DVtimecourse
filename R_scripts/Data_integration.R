## For test
#!/usr/bin/env Rscript
# Download Library
library(Seurat)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(parallelDist)
library(dplyr)
library(ComplexHeatmap)
library(gprofiler2)
options(future.globals.maxSize= 12800000000)

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

# Run standard pre-processing on Seurat object
Individual_sample_list <- readRDS(file = "/home/icbs_shared_storage_yod/Jantarika/10X_PBMC_29072018/Downstream_Analysis/Analysis_Pipeline_8/Obj/Integration/Input_SeuratObjList_10Samples_After_AddingMetadata/Input_SeuratObjList_10Samples_After_AddingMetadata_Spliced.rds")
for (RN in 1:length(Individual_sample_list)) {
  print(paste("Running sample" , RN , sep = " "))
  Individual_sample_list[[RN]] <- Standard_PreProcessing_QC(Individual_sample_list[[RN]])
}

# Integration processes
Integration_features <- SelectIntegrationFeatures(object.list = Individual_sample_list, nfeatures = 3000)
Individual_sample_list <- PrepSCTIntegration(object.list = Individual_sample_list, anchor.features = Integration_features, verbose = FALSE)
anchors <- FindIntegrationAnchors(object.list = Individual_sample_list, normalization.method = "SCT", anchor.features = Integration_features, verbose = FALSE)
sc_integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = FALSE)
sc_integrated <- RunPCA(sc_integrated, verbose = FALSE)
sc_integrated <- RunUMAP(sc_integrated, reduction = "pca", dims = 1:30)
sc_integrated <- FindNeighbors(sc_integrated, dims = 1:30)
sc_integrated <- FindClusters(sc_integrated, resolution = 3  , algorithm = 2)

# Normalize integrated data
DefaultAssay(sc_integrated) <- "RNA"
sc_integrated <- NormalizeData(sc_integrated, verbose = FALSE)

# Save integration objects
saveRDS(sc_integrated , file = "/home/icbs_shared_storage_yod/Jantarika/10X_PBMC_29072018/Downstream_Analysis/Analysis_Pipeline_8/Obj/test_script_Github/Integration/sc_integrated.rds")

