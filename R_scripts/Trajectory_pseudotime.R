#!/usr/bin/env Rscript
# Download libraries
library(tidyverse)
library(Seurat)
library(monocle3)


### READ ME ####
# This script is for running trajectories and pseudotime analyses using Monocle3
# Integrated data of each subpopulation is used as input
# Note that : Healthy PBMC samples are not included in the integration step

# Import integrated 
sc_integrated <- readRDS(file = "/PATH_TO_DIRECTORY/OBJECT_NAME.rds")

# Extract expression matrix, cell metadata and gene names from Seurat's obj
expression_matrix <- Seurat::GetAssayData(sc_integrated, assay = "RNA", slot = "data")
meta_data <- sc_integrated@meta.data
seurat_genes <- data.frame(gene_short_name = rownames(sc_integrated[["RNA"]]), row.names = rownames(sc_integrated[["RNA"]]))

# Create Monocle's Obj.
new_cds <- monocle3::new_cell_data_set(expression_data = expression_matrix, cell_metadata = meta_data, gene_metadata = seurat_genes)

# Creat SingleCellExperiment (sce) object fom Seurat obj.
sce <- as.SingleCellExperiment(sc_integrated, assay = "RNA")

# Add all Seurat reductions (PCA, UMAP) into monocle's Obj. (cds)
SingleCellExperiment::reducedDims(new_cds) <- SingleCellExperiment::reducedDims(sce)

# Import counts and data from Seurat Obj. into cds 
SummarizedExperiment::assays(new_cds) <- SummarizedExperiment::assays(sce)

# Load gene names from in Seurat into CDS
SummarizedExperiment::rowData(new_cds) <- SummarizedExperiment::rowData(sce)
SummarizedExperiment::rowData(new_cds)$gene_short_name <-  row.names(new_cds)

# Loading in Seurat gene loadings into CDS
new_cds@preprocess_aux$gene_loadings <- sc_integrated@reductions[["pca"]]@feature.loadings

# The next two commands are run in order to allow "order_cells" to be run in monocle3
rownames(new_cds@principal_graph_aux[['UMAP']]$dp_mst) <- NULL
colnames(SingleCellExperiment::reducedDims(new_cds)[["UMAP"]]) <- NULL

# Run trajectory and pseudotime
new_cds = cluster_cells(new_cds,reduction_method = "UMAP",cluster_method = "louvain")
new_cds <- learn_graph(new_cds)
new_cds <- order_cells(new_cds)

# Visualize pseudotime on UMAP
plot_cells(new_cds,
           color_cells_by = "pseudotime",
           label_cell_groups=F,
           label_leaves=F,
           label_branch_points=F,
           graph_label_size=3,
           label_roots = F,
           label_groups_by_cluster = F)


## Visualize the expression of interest gene across pseudotime
# - Extract pseudotime values
pT<-data.frame(pseudotime_coordinate=pseudotime(new_cds))

# - Extract cell types, times, expression matrices 
cell_types  <- sc_integrated@meta.data %>% data.frame %>% rownames_to_column() %>% select( rowname , Cell_Types )
times <- sc_integrated@meta.data %>% data.frame %>% rownames_to_column() %>% select( rowname , time )
expression_table <-GetAssayData(object = sc_integrated, assay = "RNA" , slot = "data") %>% data.frame %>% rownames_to_column()

# - Select genes of interest
selected_gene_table <- expression_table%>%filter(rowname %in% c("GLG1", "CCR10", "SELPLG", "CCR9", "ITGB7" , "CCR7", "CXCR4", "SELL" , "ITGAE", "ITGB2" , "CCR4" , "CCR6" , "CCR10"  , "CXCR3" , "CCR2")) %>% column_to_rownames() %>% t %>% data.frame() %>%rownames_to_column

expression_table_with_pseudotime <-pT%>%rownames_to_column%>%left_join(selected_gene_table,by="rowname")%>%left_join(cell_types,by="rowname") %>%left_join(times,by="rowname")

# Plot
ggplot(expression_table_with_pseudotime, aes(x=pseudotime_coordinate,y=GLG1  , color=Cell_Types))  + geom_point(size=1) + geom_smooth(method='loess', color="black") + xlab("Pseudotime") + ylab("Expression value") + theme_classic()   + theme(plot.title = element_text( size=15, face="bold") ,axis.text = element_text(size = 12 , colour = "black" )  , axis.title.y = element_text(color="black", size=12) , axis.title.x = element_text(colour = "black" , size = 12),legend.title = element_text(face="bold" ,size = 15), legend.text = element_text(size = 15 )) +  labs(color = "Cell types")












