#!/usr/bin/env Rscript
# Download libraries
library(Seurat)
library(SoupX)

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
  input <- FindClusters(input,  verbose = FALSE)
  
  # Run non-linear dimensional reduction (UMAP/tSNE)
  print("Run non-linear dimensional reduction (UMAP/tSNE)")
  input <- RunUMAP(input, dims = 1:30 ,  verbose = FALSE)
  return(input)
}


# Create SoupChannel objects 
# Load 10x data (raw_feature_bc_matrix and filtered_feature_bc_matrix) processed by cellRanger 
SoupXObj_list <-list()
sample_label <- c( "T002D1" , "T002D2"  , "T002D0" , "T002F1" , "I076D2" , "I076D3" ,"I076D0" , "I076F1" )
for (RN in sample_label) {
  print(RN)
  dataDirs <- paste0("/home/icbs_shared_storage_yod/Jantarika/10X_PBMC_29072018/Mapping_V3_Output/Latest_15Nov2019/" , RN ,"/" ,  RN , "/outs")
  SoupXObj_list[[RN]] <- load10X(dataDirs)
}

# Create Seurat objects 
# Load 10x data (filtered_feature_bc_matrix) processed by cellRanger
Seurat_Obj_list <- list()
sample_label <- c( "T002D1" , "T002D2"  , "T002D0" , "T002F1" , "I076D2" , "I076D3" ,"I076D0" , "I076F1" )
for (RN in sample_label) {
  print(RN)
  data_dir <- paste0("/home/icbs_shared_storage_yod/Jantarika/10X_PBMC_29072018/Mapping_V3_Output/Latest_15Nov2019/" , RN ,"/" ,  RN , "/outs/filtered_feature_bc_matrix/")
  list.files(data_dir)
  temp <- Read10X(data.dir = data_dir)
  Seurat_Obj_list[[RN]] <- CreateSeuratObject(counts = temp)
}

# Run standard pre-processing in Seurat objects 
for (RN in 1:length(Seurat_Obj_list)) {
  print(RN)
  input <- Seurat_Obj_list[[RN]]
  temp <- Standard_PreProcessing_QC(input)
  Seurat_Obj_list[[RN]] <- temp
}

# Extract clustering information from Seurat objects
DR_list <- list()
sample_label <- c( "T002D1" , "T002D2"  , "T002D0" , "T002F1" , "I076D2" , "I076D3" ,"I076D0" , "I076F1" )
for (RN in 1:length(Seurat_Obj_list)) {
  temp  <- as.data.frame(Embeddings(object = Seurat_Obj_list[[RN]][["umap"]]))
  colnames(temp ) = c('RD1','RD2')
  temp$Cluster <- factor(Seurat_Obj_list[[RN]]@meta.data$SCT_snn_res.0.8)
  DR_list[[sample_label[RN]]] <- temp
}

## Provid list of Immunoglobulin (Ig) genes
igGenes <-  c("IGLC2","IGLC3","IGKC","IGHG1","IGHG4", "JCHAIN" , "IGHA1", "IGHM" , "IGLL5" , "IGHG3")

# Estimating Non-expressing Cells 
uesToEst_list <- list()
for (RN in 1:length(SoupXObj_list)) {
  uesToEst_list[[RN]] <- estimateNonExpressingCells(SoupXObj_list[[RN]], nonExpressedGeneList = list(IG = igGenes), 
                                                    clusters = setNames(DR_list[[RN]]$Cluster, rownames(DR_list[[RN]])))
}

# Calculate the contamination fraction
for (RN in 1:length(SoupXObj_list)) {
  SoupXObj_list[[RN]] <- calculateContaminationFraction(SoupXObj_list[[RN]], list(IG = igGenes), useToEst = uesToEst_list[[RN]])
}

# Correct the expression profile 
Corrected_SoupX_list <- list()
sample_label <- c( "T002D1" , "T002D2"  , "T002D0" , "T002F1" , "I076D2" , "I076D3" ,"I076D0" , "I076F1" )
for (RN in 1:length(SoupXObj_list)) {
  print(RN)
  Corrected_SoupX_list[[RN]] <- adjustCounts(SoupXObj_list[[RN]])
}

# Save output files 
DropletUtils:::write10xCounts("/home/icbs_shared_storage_yod/Jantarika/10X_PBMC_29072018/Downstream_Analysis/Analysis_Pipeline_8/Obj/SoupX/Output_files/T002D1", Corrected_SoupX_list[[1]])
DropletUtils:::write10xCounts("/home/icbs_shared_storage_yod/Jantarika/10X_PBMC_29072018/Downstream_Analysis/Analysis_Pipeline_8/Obj/SoupX/Output_files/T002D2", Corrected_SoupX_list[[2]])
DropletUtils:::write10xCounts("/home/icbs_shared_storage_yod/Jantarika/10X_PBMC_29072018/Downstream_Analysis/Analysis_Pipeline_8/Obj/SoupX/Output_files/T002D0", Corrected_SoupX_list[[3]])
DropletUtils:::write10xCounts("/home/icbs_shared_storage_yod/Jantarika/10X_PBMC_29072018/Downstream_Analysis/Analysis_Pipeline_8/Obj/SoupX/Output_files/T002F1", Corrected_SoupX_list[[4]])
DropletUtils:::write10xCounts("/home/icbs_shared_storage_yod/Jantarika/10X_PBMC_29072018/Downstream_Analysis/Analysis_Pipeline_8/Obj/SoupX/Output_files/I076D2", Corrected_SoupX_list[[5]])
DropletUtils:::write10xCounts("/home/icbs_shared_storage_yod/Jantarika/10X_PBMC_29072018/Downstream_Analysis/Analysis_Pipeline_8/Obj/SoupX/Output_files/I076D3", Corrected_SoupX_list[[6]])
DropletUtils:::write10xCounts("/home/icbs_shared_storage_yod/Jantarika/10X_PBMC_29072018/Downstream_Analysis/Analysis_Pipeline_8/Obj/SoupX/Output_files/I076D0", Corrected_SoupX_list[[7]])
DropletUtils:::write10xCounts("/home/icbs_shared_storage_yod/Jantarika/10X_PBMC_29072018/Downstream_Analysis/Analysis_Pipeline_8/Obj/SoupX/Output_files/I076F1", Corrected_SoupX_list[[8]])

# Save SoupChannel objects
saveRDS(SoupXObj_list , file = "/home/icbs_shared_storage_yod/Jantarika/10X_PBMC_29072018/Downstream_Analysis/Analysis_Pipeline_8/Obj/SoupX/SoupXObj_list.rds")
saveRDS(Corrected_SoupX_list , file = "/home/icbs_shared_storage_yod/Jantarika/10X_PBMC_29072018/Downstream_Analysis/Analysis_Pipeline_8/Obj/SoupX/Corrected_SoupX_list.rds")

# Save Seurat objects
saveRDS(Seurat_Obj_list , file = "/home/icbs_shared_storage_yod/Jantarika/10X_PBMC_29072018/Downstream_Analysis/Analysis_Pipeline_8/Obj/SoupX/Seurat_Obj_list_PreProcessQC.rds")




