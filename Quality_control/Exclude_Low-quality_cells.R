#!/usr/bin/env Rscript
# Download libraries
library(Seurat)
library(DoubletFinder)

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


# Create Seurat objects 
# Load input genarated by SoupX
Seurat_Obj_list <- list()
sample_label <- c( "T002D1" , "T002D2"  , "T002D0" , "T002F1" , "I076D2" , "I076D3" ,"I076D0" , "I076F1" , "Healthy_3P_Our" , "Healthy_3P_10X"  )
for (RN in sample_label) {
  print(RN)
  data_dir <- paste0("/home/icbs_shared_storage_yod/Jantarika/10X_PBMC_29072018/Downstream_Analysis/Analysis_Pipeline_8/Obj/SoupX/Output_files/" , RN )
  print(data_dir)
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

# 1.Exclude cells that have more than 10% of mitochondrial counts
Seurat_Obj_AfterQC_CutPercentMT_list <- list()
for (RN in 1:length(Seurat_Obj_list)) {
  Seurat_Obj_AfterQC_CutPercentMT_list[[RN]] <- subset(Seurat_Obj_list[[RN]], subset =  percent.mt < 10)
}


# Create Seurat objects 
# Load input genarated by SoupX
Seurat_Obj_list <- list()
sample_label <- c( "T002D1" , "T002D2"  , "T002D0" , "T002F1" , "I076D2" , "I076D3" ,"I076D0" , "I076F1" , "Healthy_3P_Our" , "Healthy_3P_10X"  )
for (RN in sample_label) {
  print(RN)
  data_dir <- paste0("/home/icbs_shared_storage_yod/Jantarika/10X_PBMC_29072018/Downstream_Analysis/Analysis_Pipeline_8/Obj/SoupX/Output_files/" , RN )
  print(data_dir)
  list.files(data_dir)
  temp <- Read10X(data.dir = data_dir)
  Seurat_Obj_list[[RN]] <- CreateSeuratObject(counts = temp)
}


# Extract only cells that have more than 10% of mitochondrial counts  
Seurat_Obj_after_Cutting10PMito <- list()
for (RN in 1:length(Seurat_Obj_list)) {
  Seurat_Obj_after_Cutting10PMito[[RN]] <- subset(Seurat_Obj_list[[RN]] , cells = colnames(Seurat_Obj_AfterQC_CutPercentMT_list[[RN]]))
}
names(Seurat_Obj_after_Cutting10PMito) <- names(Seurat_Obj_list)

# Run standard pre-processing in Seurat objects 
for (RN in 1:length(Seurat_Obj_after_Cutting10PMito)) {
  input <- Seurat_Obj_after_Cutting10PMito[[RN]]
  temp <- Standard_PreProcessing_QC(input)
  Seurat_Obj_after_Cutting10PMito[[RN]] <- temp
}


# 2. Run DoubletFinder
for (i in 1:length(Seurat_Obj_after_Cutting10PMito)) {
  print(paste0("start ",i))
  print("pK Identification (no ground-truth)")
  sweep.res.list_temp <- paramSweep_v3(Seurat_Obj_after_Cutting10PMito[[i]], PCs = 1:30, sct = TRUE,num.cores = 25)
  sweep.stats_temp <- summarizeSweep(sweep.res.list_temp, GT = FALSE)
  pK_object <- find.pK(sweep.stats_temp)
  pK_value<-pK_object%>%filter(MeanBC==max(MeanBC))%>%dplyr::select(pK)%>%unlist%>%as.character()%>%as.numeric()
  
  print("Homotypic Doublet Proportion Estimate")
  annotations <- Seurat_Obj_after_Cutting10PMito[[i]]@meta.data$seurat_clusters
  annotations <- as.character(annotations)
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(0.039*length(Seurat_Obj_after_Cutting10PMito[[i]]%>%colnames))  ## Assuming 7.5% doublet formation rate - tailor for your dataset!!##!!
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  ##pK<0.1
  print("Run DoubletFinder with varying classification stringencies")
  Seurat_Obj_after_Cutting10PMito[[i]]  <- doubletFinder_v3(Seurat_Obj_after_Cutting10PMito[[i]] , PCs = 1:30, pN = 0.25, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
  print(paste0("finish ",i))
}

# Save Seurat object without cutting doublet 
saveRDS(Seurat_Obj_after_Cutting10PMito , file = "/home/icbs_shared_storage_yod/Jantarika/10X_PBMC_29072018/Downstream_Analysis/Analysis_Pipeline_8/Obj/DoubletFinder/Output/Seurat_Obj_SoupX_Cut10PercentMT_DoubletCalculation_list.rds")






