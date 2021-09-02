# scRNAseq_DENVtimecourse Project
# Part 02 : Data normalisation and integration
# Author : Jantarika Kumar Arora


# ------------------------------------------
# README
# ------------------------------------------
# This script represents the steps of data normalisation and integration using Seurat package v.3.1.2
# The processed data generated from part 01 were used as the input


# ------------------------------------------
# Load required libraries 
# ------------------------------------------
library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)

options(future.globals.maxSize= 12800000000)

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
  input <- FindClusters(input,  verbose = FALSE , algorithm = 2)
  
  # Run non-linear dimensional reduction (UMAP)
  print("Run non-linear dimensional reduction (UMAP)")
  input <- RunUMAP(input, dims = 1:30 ,  verbose = FALSE)
  return(input)
}

# Import processed data generated from part01
Individual_sample_list <- readRDS(file = "/PATH_TO_DIRECTORY/OBJECT_NAME.rds") # Directory of processed data generated from part01

# Run standard pre-processing on Seurat object
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

# Normalise integrated data on RNA assay
DefaultAssay(sc_integrated) <- "RNA"
sc_integrated <- NormalizeData(sc_integrated, verbose = FALSE)

# Save integration object
saveRDS(sc_integrated , file = "/PATH_TO_DIRECTORY/sc_integrated.rds")

