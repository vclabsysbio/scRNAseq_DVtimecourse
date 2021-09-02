# scRNAseq_DENVtimecourse Project
# Part03 : Highly variable genes (HVGs) and biological process (BP) analyses (Figure 2A-C)
# Author : Jantarika Kumar Arora


# ------------------------------------------
# README
# ------------------------------------------
# This script shows how to calculate the pearson correlation, HVG, heatmap, and BP analyses.
# The integrated data generated from part02 was used as the input


# ------------------------------------------
# Load required libraries 
# ------------------------------------------
library(Seurat)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(parallelDist)
library(dplyr)
library(ComplexHeatmap)
library(gprofiler2)


# Import integrated data
sc_integrated <- readRDS(file = "/PATH_TO_DIRECTORY/sc_integrated.rds") # Directory of processed data generated from part02


# Calculate the average gene expression of all cells
Idents(sc_integrated) <- "severity_time" # can be changed by user
Average_data <- AverageExpression(sc_integrated , assays = "RNA" , slot = "data" )

# Perform Principal Component Analysis (PCA) (Figure 2A)
pca <- prcomp(t(Average_data$RNA))
pca_perc <- round(100*pca$sdev^2/sum(pca$sdev^2),1)
df_pca <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], sample = colnames(Average_data$RNA))
# - Add metadata
df_pca$Severity <- c("DF" , "DF" , "DF" , "DF" , "DHF" , "DHF" , "DHF", "DHF"  , "Healthy" , "Healthy") # Can be changed by user
df_pca$Time <- c("Day-2" , "Day-1" , "Def" , "Wk2" , "Day-2" , "Day-1" , "Def" , "Wk2"  ,"Healthy I" , "Healthy II") # Can be changed by user
df_pca$Time <- factor(df_pca$Time , levels = c( "Day-2" , "Day-1" , "Def" , "Wk2" ,  "Healthy I" , "Healthy II")) # Can be changed by user

# Plot PCA   
pca_plot <- ggplot(df_pca, aes(PC1,PC2, color = Time))+ geom_point(aes(shape = Severity ), size=6 , stroke = 1.4)+ labs(x=paste0("PC1 (",pca_perc[1],"%)"), y=paste0("PC2 (",pca_perc[2],"%)"))  + theme(axis.text = element_text(size = 17 , face="bold" , colour = "black") ,  axis.title.y = element_text(color="black", size=15, face="bold") , axis.title.x  = element_text(color="black", size=17, face="bold") , legend.title = element_text(face = "bold" , size = 17)  ,  legend.text = element_text(size = 16) , legend.key.size = unit(1, "cm") , legend.key.width = unit(0.5,"cm") , legend.key = element_rect(fill = "white") ) + scale_shape_manual(values = c(7,16 , 2 , 17)) 
plot(pca_plot)

# Heatmap (Figure 2B)
# Extract average expression of 500 top genef from PC1 and PC2 (HVGs)
PC1_Genes <- data.frame(sort(abs(pca$rotation[,"PC1"]), decreasing=TRUE)[1:500])
PC2_Genes <- data.frame(sort(abs(pca$rotation[,"PC2"]), decreasing=TRUE)[1:500])
union_genes <- union(rownames(PC1_Genes) , rownames(PC2_Genes))
input_heatmap <- Average_data$RNA[rownames(Average_data$RNA) %in% union_genes,]

# Normalise by Z-score
input_heatmap <- t(apply((input_heatmap), 1, function(x){
  mean <- mean(x)
  SD <- sd(x)
  Z_score <- (x-mean)/SD
  Z_score
}))

# Construct heatmap
input_heatmap <- as.data.frame(input_heatmap)
parallelDist_dtw <- hclust(parDist(x = as.matrix(input_heatmap), method = "dtw") , method="ward")
parallelDist_dtw_log2_de_or <-input_heatmap[parallelDist_dtw$order,]

# Elbow plot to determine the number of clusters
ggplot(parallelDist_dtw$height %>% as.tibble() %>% add_column(groups = length(parallelDist_dtw$height):1) , aes(x=groups, y=value)) +geom_point() + geom_line() + xlim(0,20)

row_cutree <- data.frame(cutree(parallelDist_dtw, k = 8 )) 
colnames(row_cutree) <- "cluster"
parallelDist_dtw_log2_de_or <- parallelDist_dtw_log2_de_or %>%  rownames_to_column %>%  right_join( row_cutree %>%  rownames_to_column %>%  select(cluster ,rowname ) , by = "rowname")
parallelDist_dtw_log2_de_or <- parallelDist_dtw_log2_de_or[parallelDist_dtw$order,]%>% cbind.data.frame(index=seq(1:nrow(.)))
DF_cluster_annotation <- rowAnnotation( df = parallelDist_dtw_log2_de_or%>%remove_rownames()%>%column_to_rownames()%>%select(cluster) ,  name = "cluster" , show_annotation_name = c(TRUE) , width = unit (0.5 , "cm") , show_legend = TRUE , col = list(cluster = c( "1" = "purple4" , "2" = "royalblue1"  , "3" = "orangered2" , "4" = "chartreuse4" ,  "5" = "hotpink2" , "6" = "mediumaquamarine" , "7" = "midnightblue" , "8" = "lightslateblue"  )))
out_heatmap <- Heatmap(parallelDist_dtw_log2_de_or%>%remove_rownames()%>%column_to_rownames()%>%select(-cluster,-index) , col = circlize::colorRamp2(c(-3,0,3), c("green","black","red")), name = "Average Gene Expression (Z-score)",  cluster_columns = F, show_column_dend = F, cluster_rows = F, show_row_names = F, show_column_names = F, na_col = "grey" , row_names_gp = gpar(fontsize = 1))

# Visualise
draw(out_heatmap+DF_cluster_annotation, row_split = parallelDist_dtw_log2_de_or$cluster)


# BP analyses
# Create the input
# Note that : the number of clusters is corresponding to numbers of clusters on heatmap
input_GO <-list(cluster1=row_cutree %>% rownames_to_column() %>% filter(cluster=="1")%>%select(rowname)%>%unlist%>%as.character(), cluster2=row_cutree %>% rownames_to_column() %>% filter(cluster=="2")%>%select(rowname)%>%unlist%>%as.character(), cluster3=row_cutree %>% rownames_to_column() %>% filter(cluster=="3")%>%select(rowname)%>%unlist%>%as.character(), cluster4=row_cutree %>% rownames_to_column() %>% filter(cluster=="4")%>%select(rowname)%>%unlist%>%as.character(), cluster5=row_cutree %>% rownames_to_column() %>% filter(cluster=="5")%>%select(rowname)%>%unlist%>%as.character(),
                cluster6=row_cutree %>% rownames_to_column() %>% filter(cluster=="6")%>%select(rowname)%>%unlist%>%as.character(),                   
                cluster7=row_cutree %>% rownames_to_column() %>% filter(cluster=="7")%>%select(rowname)%>%unlist%>%as.character(),
                cluster8=row_cutree %>% rownames_to_column() %>% filter(cluster=="8")%>%select(rowname)%>%unlist%>%as.character())
bg <- rownames(Average_data$RNA)

out_GO <- list()
for (RN in 1:length(input_GO)) {
  out_GO[[paste0("Cluster" , RN)]] <- gost(query = input_GO[[RN]] , organism = "hsapiens" , correction_method = "fdr" , custom_bg = rownames(bg) , significant = TRUE , user_threshold = 0.05 , evcodes = TRUE)
}

