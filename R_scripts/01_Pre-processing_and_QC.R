# scRNAseq_DENVtimecourse Project
# Part 01 : Pre-processing steps and quality controls
# author : Jantarika Kumar Arora


# ------------------------------------------
# README
# ------------------------------------------
# This script covers the pre-processing steps and quality controls in each sample. 
# The raw and filtered expression matrices from cellranger count were used as the input.
# Step 1 : Correct the expression matrices by SoupX 
# Step 2 : Discard cells that have more than 10% of mitochondrial counts
# Step 3 : Exclude doublets by DoubletFinder



# ------------------------------------------
# Load required libraries 
# ------------------------------------------
library(Seurat)
library(SoupX)
library(DropletUtils)
library(DoubletFinder)
library(stringr)
library(dplyr)
library(tidyverse)

# ------------------------------------------
# Step 1 : Correct the expression matrices by SoupX 
# ------------------------------------------

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
sc_SoupX <- load10X("/PATH_TO_DIRECTORY/") # Directory of raw_ and filtered_feature_bc_matrix processed by cellRangercount

# Create Seurat object
sc_Seurat <- Read10X(data.dir = "/PATH_TO_DIRECTORY/") # Directory of filtered_feature_bc_matrix processed by cellRangercount
sc_Seurat <- CreateSeuratObject(counts = sc_Seurat)

# Run standard pre-processing on Seurat object
sc_Seurat <- Standard_PreProcessing_QC(sc_Seurat)

# Add clustering information 
clustering_info <- data.frame(Cluster = sc_Seurat@meta.data$SCT_snn_res.0.8)
rownames(clustering_info) <- rownames(sc_Seurat@meta.data)
sc_SoupX <-  setClusters(sc_SoupX, setNames(clustering_info$Cluster, rownames(clustering_info)))

# Provide the list of Immunoglobulin (Ig) genes.
igGenes <-  c("IGLC2","IGLC3","IGKC","IGHG1","IGHG4", "JCHAIN" , "IGHA1", "IGHM" , "IGLL5" , "IGHG3")

# Estimate non-expressing cells 
useToEst <- estimateNonExpressingCells(sc_SoupX, nonExpressedGeneList = list(IG = igGenes))

# Calculate the contamination fraction
sc_SoupX <- calculateContaminationFraction(sc_SoupX, list(IG = igGenes), useToEst = useToEst)

# Correct the expression profile 
sc_SoupX_output <- adjustCounts(sc_SoupX)

# Save corrected expression matrices 
DropletUtils:::write10xCounts("/PATH_TO_DIRECTORY/", sc_SoupX_output) 

# ------------------------------------------
# Step 2 : Discard cells that have more than 10% of mitochondrial counts
# ------------------------------------------

# Create Seurat object
sc_Seurat <- Read10X(data.dir = "/PATH_TO_DIRECTORY/") # Directory of corrected expression matrices processed by SoupX
sc_Seurat <- CreateSeuratObject(counts = sc_Seurat)
sc_Seurat_procesing <- sc_Seurat

# Run standard pre-processing on Seurat object
sc_Seurat_procesing <- Standard_PreProcessing_QC(sc_Seurat_procesing)

# Exclude cells that have more than 10% of mitochondrial counts
sc_Seurat_procesing <- subset(sc_Seurat_procesing, subset =  percent.mt < 10)

# Extract the cells that pass % mitochondrial cut-off
sc_Seurat_procesing <- subset(sc_Seurat , cells = colnames(sc_Seurat_procesing))

# ------------------------------------------
# Step 3 : Exclude doublets by DoubletFinder
# ------------------------------------------

# Run standard pre-processing on Seurat object with cells passeing % mitochondrial cut-off 
sc_Seurat_procesing <- Standard_PreProcessing_QC(sc_Seurat_procesing)

# Run DoubletFinder
sweep.res.list <- paramSweep_v3(sc_Seurat_procesing, PCs = 1:30, sct = TRUE,num.cores = 25)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
pK_value <- find.pK(sweep.stats)

# Extract the pK value at maximum mean-variance-normalized bimodality coefficient (MeanBC)  
pK_value<-pK_value%>%filter(MeanBC==max(MeanBC))%>%dplyr::select(pK)%>%unlist%>%as.character()%>%as.numeric()

# Predict doublets 
annotations <- sc_Seurat_procesing@meta.data$seurat_clusters
annotations <- as.character(annotations)
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.039*length(sc_Seurat_procesing%>%colnames))  # Assuming 3.9% doublet formation rate according to the 10x protocol
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
sc_Seurat_procesing <- doubletFinder_v3(sc_Seurat_procesing , PCs = 1:30, pN = 0.25, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
colnames(sc_Seurat_procesing@meta.data)[length(sc_Seurat_procesing@meta.data)] <- "DoubletFinder"

# Visualise the doublets 
DimPlot(sc_Seurat_procesing , group.by = "DoubletFinder" , cols = c("red" , "black"))

# Exclude the doublets
sc_Seurat <- subset(sc_Seurat , cells =  rownames(sc_Seurat_procesing@meta.data[sc_Seurat_procesing@meta.data$DoubletFinder == "Singlet",])) 

# Add necessary information into metadata for data integration in the next step
sc_Seurat$time <- "time" # Time-course of DENV infection eg "Day-2"
sc_Seurat$severity <- "severity" # Severity of patient eg "DF"
sc_Seurat$severity_time <- paste(sc_Seurat$severity, sc_Seurat$severity)

# Save object for the next step, data integration. 
saveRDS(sc_Seurat , file = "/PATH_TO_DIRECTORY/OBJECT_NAME.rds")


