#!/usr/bin/env Rscript
# Download Library
library(Seurat)

# Import individual samples with cells that pass qc
Input_SeuratObjList_10Samples_After_AddingMetadata_Spliced <- readRDS(file = "/home/icbs_shared_storage_yod/Jantarika/10X_PBMC_29072018/Downstream_Analysis/Analysis_Pipeline_8/Obj/Integration/Input_SeuratObjList_10Samples_After_AddingMetadata/Input_SeuratObjList_10Samples_After_AddingMetadata_Spliced.rds")

# Standard pre-processing of Seurat object including SCT normalization, PCs = 1:30
Standard_PreProcessing_QC <- function(input){
  # cal % mito
  print("cal % mito")
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
  input <- FindClusters(input,  verbose = FALSE , algorithm = 2 )
  
  # Run non-linear dimensional reduction (UMAP/tSNE)
  print("Run non-linear dimensional reduction (UMAP/tSNE)")
  input <- RunUMAP(input, dims = 1:30 ,  verbose = FALSE )
  return(input)
}

# Run standard pre-processing in Seurat objects
for (RN in 1:length(Input_SeuratObjList_10Samples_After_AddingMetadata_Spliced)) {
  print(RN)
  input <- Input_SeuratObjList_10Samples_After_AddingMetadata_Spliced[[RN]]
  temp <- Standard_PreProcessing_QC(input)
  Input_SeuratObjList_10Samples_After_AddingMetadata_Spliced[[RN]] <- temp
}

# Integration process
Integration_features <- SelectIntegrationFeatures(object.list = Input_SeuratObjList_10Samples_After_AddingMetadata_Spliced, nfeatures = 3000)
Input_SeuratObjList_10Samples_After_AddingMetadata_Spliced <- PrepSCTIntegration(object.list = Input_SeuratObjList_10Samples_After_AddingMetadata_Spliced, anchor.features = Integration_features, verbose = FALSE)

# Find anchors
anchors <- FindIntegrationAnchors(object.list = Input_SeuratObjList_10Samples_After_AddingMetadata_Spliced, normalization.method = "SCT", 
                                  anchor.features = Integration_features, verbose = FALSE)
# Integrate Data
Integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = FALSE)
Integrated <- RunPCA(Integrated, verbose = FALSE)
Integrated <- RunUMAP(Integrated, reduction = "pca", dims = 1:30)
Integrated <- FindNeighbors(Integrated, dims = 1:30)
Integrated <- FindClusters(Integrated, resolution = 1.5  , algorithm = 2)

# Save Integrated object
saveRDS(Integrated , file = "/home/icbs_shared_storage_yod/Jantarika/10X_PBMC_29072018/Downstream_Analysis/Analysis_Pipeline_8/Obj/Integration/FindCluster_LouvainMultiRefinement/Integrated10Samples_LouvainMultiRefinement_for_FindCluster_SCTransform.rds")

