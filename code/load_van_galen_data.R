library(exFINDER)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(RColorBrewer)
library(dplyr)
library(ComplexHeatmap)
library(ggraph)
library(stringr)
library(igraph)

data = load("AML1 data.RData")

## Plot UMAP based on log normalized counts

ElbowPlot(AML) #to determine dimensionality
AML <- FindNeighbors(AML, reduction = "pca", dims = 1:20)
AML <- FindClusters(AML,verbose = FALSE, resolution=0.2) 
AML <- RunUMAP(AML, dims = 1:20)
DimPlot(AML, reduction = "umap",  group.by = 'cell_type', label= TRUE)
DimPlot(AML, reduction = "umap",  group.by = c("cell_type", 'cell_subtype'), label= TRUE)

## Check that the counts are normalized

AML@assays$SCT$scale.data %>% as_tibble()
AML@assays$SCT$counts %>% str
AML@assays$SCT$counts[1:100, 1:100] %>% as.matrix() %>% as.vector() %>% unique

## Plots UMAP based on raw counts
# # AML <- SCTransform(AML, verbose = TRUE) # run sctransform
# AML <- FindNeighbors(AML, reduction = "pca", dims = 1:20)
# AML <- RunUMAP(AML, reduction = "pca", dims = 1:20)


meta_data = AML@meta.data %>% as_tibble()

meta_data$cell_type %>% unique
meta_data$cell_subtype %>% unique %>% length
meta_data$cell_cycle_phase %>% unique

ggplot(data = meta_data, aes(umap1, umap2, color = cell_type)) + geom_point()
ggplot(data = meta_data, aes(umap1, umap2, color = cell_subtype)) + geom_point()

table(meta_data$cell_type, meta_data$cell_subtype) %>% lattice::levelplot()

## select character variables in meta_data
meta_data %>% select_if(is.character) %>% map(unique)

table(meta_data$source, meta_data$cell_subtype)


## Check if the data contains the 6 healty samples -----
## Read csv file
cell_info_df = read.csv("/Users/vgoepp/work/aachen/laila/LailaProject/data/Cells.csv", row.names=1) %>% 
  as_tibble()

## For each variable in cell_info_df, if it is not a numeric variable, print its unique values
cell_info_df %>% select_if(is.character) %>% map(unique)


cell_line_df = cell_info_df %>% filter(source == "cell_line")

## vanGalen postprocessed data -----
library(tidyverse)
vanGalen = readRDS("/Users/vgoepp/work/aachen/laila/LailaProject/data_whole_processed/Seurat_AML.rds")
vanGalen@assays$RNA$data %>% str
vanGalen@assays$RNA$counts %>% str

## Convert counts to an integer-valued sparse matrix
sparse_counts = vanGalen@assays$RNA$counts %>% as("dgCMatrix")
print(object.size(sparse_counts) / object.size(vanGalen@assays$RNA$counts))
vanGalen@assays$RNA$counts = sparse_counts

sparse_data = vanGalen@assays$RNA$data %>% as("dgCMatrix")
print(object.size(sparse_data) / object.size(vanGalen@assays$RNA$data))
vanGalen@assays$RNA$data = sparse_data

## Save the data with sparse assays
saveRDS(vanGalen, "/Users/vgoepp/work/aachen/laila/LailaProject/data_whole_processed/Seurat_AML_sparse.rds")

vanGalen_sparse = readRDS("/Users/vgoepp/work/aachen/laila/LailaProject/data_whole_processed/Seurat_AML_sparse.rds")
object.size(vanGalen_sparse)

