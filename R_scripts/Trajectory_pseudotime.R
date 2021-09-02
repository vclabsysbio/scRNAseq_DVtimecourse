# scRNAseq_DENVtimecourse Project
# Part -- : Trajectory and pseudotime analyses
# Author : Jantarika Kumar Arora


# ------------------------------------------
# README
# ------------------------------------------
# This script represents the steps performing trajectory and pseudotime analyses by Monocle3 (v.0.2.3.0)
# Integrated data of subpopulations were used as the input
# Healthy samples were not included in the integrated data
# The Effector CD8-3 and Naive CD4 T cells were set at the root of the pseudotime for the CD8 and CD4 T cell populations

# ------------------------------------------
# Load required libraries 
# ------------------------------------------
library(tidyverse)
library(Seurat)
library(monocle3)


# Import integrated 
sc_integrated_subpops <- readRDS(file = "/PATH_TO_DIRECTORY/sc_integrated_subpops.rds") # Directory of processed data generated from part01

# Extract expression matrix, cell metadata and gene names from Seurat's obj
expression_matrix <- Seurat::GetAssayData(sc_integrated_subpops, assay = "RNA", slot = "data")
meta_data <- sc_integrated_subpops@meta.data
seurat_genes <- data.frame(gene_short_name = rownames(sc_integrated_subpops[["RNA"]]), row.names = rownames(sc_integrated_subpops[["RNA"]]))

# Create Monocle's Obj.
new_cds <- monocle3::new_cell_data_set(expression_data = expression_matrix, cell_metadata = meta_data, gene_metadata = seurat_genes)

# Creat SingleCellExperiment (sce) object fom Seurat obj.
sce <- as.SingleCellExperiment(sc_integrated_subpops, assay = "RNA")

# Add all Seurat reductions (PCA, UMAP) into monocle's Obj. (cds)
SingleCellExperiment::reducedDims(new_cds) <- SingleCellExperiment::reducedDims(sce)

# Import counts and data from Seurat Obj. into cds 
SummarizedExperiment::assays(new_cds) <- SummarizedExperiment::assays(sce)

# Load gene names from in Seurat into CDS
SummarizedExperiment::rowData(new_cds) <- SummarizedExperiment::rowData(sce)
SummarizedExperiment::rowData(new_cds)$gene_short_name <-  row.names(new_cds)

# Loading in Seurat gene loadings into CDS
new_cds@preprocess_aux$gene_loadings <- sc_integrated_subpops@reductions[["pca"]]@feature.loadings

# The next two commands are run in order to allow "order_cells" to be run in monocle3
rownames(new_cds@principal_graph_aux[['UMAP']]$dp_mst) <- NULL
colnames(SingleCellExperiment::reducedDims(new_cds)[["UMAP"]]) <- NULL

# Run trajectory and pseudotime
new_cds <- cluster_cells(new_cds,reduction_method = "UMAP",cluster_method = "louvain")
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


# Visualize the expression of interest gene across pseudotime
# Extract pseudotime values
pT <-data.frame(pseudotime_coordinate=pseudotime(new_cds))

# Extract cell types, times, expression matrices 
cell_types <- sc_integrated_subpops@meta.data %>% data.frame %>% rownames_to_column() %>% select( rowname , Cell_Types )
times <- sc_integrated_subpops@meta.data %>% data.frame %>% rownames_to_column() %>% select( rowname , time )
expression_table <-GetAssayData(object = sc_integrated_subpops, assay = "RNA" , slot = "data") %>% data.frame %>% rownames_to_column()

# Select genes of interest
selected_gene_table <- expression_table%>%filter(rowname %in% c("GLG1", "CCR10", "SELPLG" , "ITGAE", "ITGB2" , "FABP2"  , "FABP5"   ,"CCR9", "ITGB7" , "CCR7", "CXCR4", "SELL" , "CXCR3" , "CCR2")) %>% column_to_rownames() %>% t %>% data.frame() %>%rownames_to_column

expression_table_with_pseudotime <- pT%>%rownames_to_column%>%left_join(selected_gene_table,by="rowname")%>%left_join(cell_types,by="rowname") %>%left_join(times,by="rowname")

# Figure 3E
color_RN <- c("orange", "red", "#EA8331" )
ggplot(expression_table_with_pseudotime, aes(x=pseudotime_coordinate,y=SELPLG  , color=Cell_Types))  + geom_point(size=1) +  ggtitle("*SELPLG* (encoding CLA)") +  geom_smooth(method='loess', color="black") + xlab("Pseudotime") + ylab("Expression value") + theme_classic() + scale_color_manual(values=color_RN)  + theme( plot.title = element_markdown(hjust = 0, vjust = 0.5, margin = ggplot2::margin(5, 0, 5, 0) , size = 30) , axis.text = element_text(size = 22 , colour = "black" )  , axis.title.y = element_text(color="black", size=27) , axis.title.x = element_text(colour = "black" , size = 27),legend.title = element_text( size = 20) )  +ylim(c(0,3.7)) + xlim(c(0,28.5))

