#Project: detecting external signals that can drive disease and fibrosis in AML
#Laila Mohamed, 2024

#load required libraries

library(Seurat)
library(Matrix)
library(ggplot2)

#Download data-------
counts<-readMM("C:\\Users\\Dell\\Downloads\\bachelor\\Exp_data_UMIcounts.mtx")
genes<- read.table("C:\\Users\\Dell\\Downloads\\bachelor\\Genes.txt", quote="\"", comment.char="")
rownames(counts)<-genes$V1
metadata <- read.csv("C:/Users/Dell/Downloads/bachelor/meta data Cells.csv", row.names=1)
colnames(counts)<-rownames(metadata)
AML <- CreateSeuratObject(counts = counts,
                          meta.data = metadata)
#visualizing umap-------
umap1 <- AML@meta.data$umap1
umap2<-AML@meta.data$umap2
umap2<-as.data.frame(umap2)
umap_df <- as.data.frame(umap1)
umap_df$umap2<-umap2$umap2
umap_df$cell_type<-AML@meta.data$cell_type
umap_df$cell_subtype<- AML@meta.data$cell_subtype

#umap based on celltype
ggplot(umap_df, aes(x = umap1, y = umap2, color= cell_type)) +
  geom_point(size=0.3) +
  theme_minimal() +
  labs(title = "Cell type AML", x = "UMAP 1", y = "UMAP 2") +
  theme(legend.position = "right")+
  guides(color = guide_legend(override.aes = list(size = 4)))

#umap based on cell subtype
ggplot(umap_df, aes(x = umap1, y = umap2, color= cell_subtype)) +
  geom_point(size=0.3) +
  theme_minimal() +
  labs(title = "Cell subtype AML", x = "UMAP 1", y = "UMAP 2") +
  theme(legend.position = "right")+
  guides(color = guide_legend(override.aes = list(size = 4)))

#QC-----
AML <- PercentageFeatureSet(AML, pattern = "^MT-", col.name = "percent.mt")# store mitochondrial percentage in object meta data
VlnPlot(AML, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) #percent mt below threshold; no filtering needed

#normalize data

options(future.globals.maxSize = 1e+09) #maximizing memory usage
AML <- SCTransform(AML, verbose = TRUE) # run sctransform
AML <- RunPCA(AML, assay = "SCT", verbose = FALSE)
ElbowPlot(AML) #to determine dimensionality
AML <- FindNeighbors(AML, reduction = "pca", dims = 1:20)
AML <- FindClusters(AML,verbose = FALSE, resolution=0.2) 
head(Idents(AML), 5)
AML <- RunUMAP(AML, dims = 1:20)
DimPlot(AML, reduction = "umap",  group.by = 'cell_type', label= TRUE)

DimPlot(AML, reduction = "umap",  group.by = 'cell_subtype', label= TRUE, label.size=3)
