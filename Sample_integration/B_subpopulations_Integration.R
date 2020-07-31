#!/usr/bin/env Rscript
# Download Library
library(Seurat)

# Import individual samples that contain only Plasma cells, Plasmablasts, and B cells
B_Subpops_Extract_from_NewGlobal_obj_Individual_sample_list_LouvainMultiLevelRefinement_for_FindCluster_SCT <- readRDS(file = "/home/icbs_shared_storage_yod/Jantarika/10X_PBMC_29072018/Downstream_Analysis/Analysis_Pipeline_8/Obj_RNANormAssay/B_Cells/New_Analysis_inputFrom_newGlobal/B_obj_input_integration/B_PlasmaCells_Plasmablasts_Subset_Individual_sample_list.rds")


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
  input <- FindClusters(input,  verbose = FALSE  , algorithm = 2 )
  
  # Run non-linear dimensional reduction (UMAP/tSNE)
  print("Run non-linear dimensional reduction (UMAP/tSNE)")
  input <- RunUMAP(input, dims = 1:30 ,  verbose = FALSE)
  return(input)
}

# Run standard pre-processing in Seurat objects
for (RN in 1:length(B_Subpops_Extract_from_NewGlobal_obj_Individual_sample_list_LouvainMultiLevelRefinement_for_FindCluster_SCT)) {
  print(RN)
  input <- B_Subpops_Extract_from_NewGlobal_obj_Individual_sample_list_LouvainMultiLevelRefinement_for_FindCluster_SCT[[RN]]
  temp <- Standard_PreProcessing_QC(input)
  B_Subpops_Extract_from_NewGlobal_obj_Individual_sample_list_LouvainMultiLevelRefinement_for_FindCluster_SCT[[RN]] <- temp
}

# Integration process
options(future.globals.maxSize= 12800000000)
Integration_features <- SelectIntegrationFeatures(object.list = B_Subpops_Extract_from_NewGlobal_obj_Individual_sample_list_LouvainMultiLevelRefinement_for_FindCluster_SCT, nfeatures = 3000)
B_Subpops_Extract_from_NewGlobal_obj_Individual_sample_list_LouvainMultiLevelRefinement_for_FindCluster_SCT <- PrepSCTIntegration(object.list = B_Subpops_Extract_from_NewGlobal_obj_Individual_sample_list_LouvainMultiLevelRefinement_for_FindCluster_SCT, anchor.features = Integration_features, verbose = FALSE)

# Find anchors
anchors <- FindIntegrationAnchors(object.list = B_Subpops_Extract_from_NewGlobal_obj_Individual_sample_list_LouvainMultiLevelRefinement_for_FindCluster_SCT, normalization.method = "SCT", 
                                  anchor.features = Integration_features, verbose = FALSE)


# Integrate Data
Integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = FALSE)
Integrated <- RunPCA(Integrated, verbose = FALSE)
Integrated <- RunUMAP(Integrated, reduction = "pca", dims = 1:30)
Integrated <- FindNeighbors(Integrated, dims = 1:30)
Integrated <- FindClusters(Integrated  , algorithm = 2 , resolution = 3)

# Save Integrated object
saveRDS(Integrated , file = "/home/icbs_shared_storage_yod/Jantarika/10X_PBMC_29072018/Downstream_Analysis/Analysis_Pipeline_8/Obj_RNANormAssay/B_Cells/New_Analysis_inputFrom_newGlobal/Integration_Output/BSubpop_Integrated10Samples_LouvainMultiLevelRefinement_for_FindCluster_SCTransform_RNANorm.rds")


