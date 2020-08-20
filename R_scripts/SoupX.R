#!/usr/bin/env Rscript
# Download libraries
library(Seurat)
library(SoupX)
library(DropletUtils)

### READ ME ####
# This script is for correcting the expression of genes by SoupX using raw and filtered expression matrices from CellRanger as input
# Raw and filtered expression matrices are provided in A-XXXX-n.   

# Note that   
# - Soup is not applied to healthy controls, there is no ambient RNA contamination (Figure S1)
# - Corrected expression matrices by SoupX are also provided in A-XXXX-n.


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


# Create SoupChannel object
sc_SoupX <- load10X("/PATH_TO_DIRECTORY/") # Directory of raw_ and filtered_feature_bc_matrix by CellRanger

# Create Seurat object
sc_Seurat <- Read10X(data.dir = "/PATH_TO_DIRECTORY/") # Directory of filtered_feature_bc_matrix by CellRanger
sc_Seurat <- CreateSeuratObject(counts = sc_Seurat)

# Run standard pre-processing on Seurat object
sc_Seurat <- Standard_PreProcessing_QC(sc_Seurat)

# Add clustering information 
clustering_info <- data.frame(Cluster = sc_Seurat@meta.data$SCT_snn_res.0.8)
rownames(clustering_info) <- rownames(sc_Seurat@meta.data)
sc_SoupX <-  setClusters(sc_SoupX, setNames(clustering_info$Cluster, rownames(clustering_info)))

# Provide the list of Immunoglobulin (Ig) genes that are used for all samples
igGenes <-  c("IGLC2","IGLC3","IGKC","IGHG1","IGHG4", "JCHAIN" , "IGHA1", "IGHM" , "IGLL5" , "IGHG3")

# Estimate non-expressing cells 
useToEst <- estimateNonExpressingCells(sc_SoupX, nonExpressedGeneList = list(IG = igGenes))

# Calculate the contamination fraction
sc_SoupX <- calculateContaminationFraction(sc_SoupX, list(IG = igGenes), useToEst = useToEst)

# Correct the expression profile 
sc_SoupX_output <- adjustCounts(sc_SoupX)

# Save corrected expression matrices 
DropletUtils:::write10xCounts("/PATH_TO_DIRECTORY/", sc_SoupX_output) 




