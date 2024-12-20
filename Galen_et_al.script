#Project: detecting external signals that can drive disease and fibrosis in AML
#Laila Mohamed, 2024

#load required libraries

library(Seurat)
library(Matrix)
library(ggplot2)
library(exFINDER)
library(patchwork)
library(RColorBrewer)
library(dplyr)
library(ComplexHeatmap)
library(ggraph)
library(stringr)
library(igraph)
library(forcats)

#Download data-------

AML<-readRDS("C:\\Users\\Dell\\Downloads\\bachelor\\Van_Galen.rds")


#checking presence of all samples and celltypes-------
unique(AML@meta.data$CellType)
unique(AML@meta.data$orig.ident)

#quality control and umaps------
AML<- subset (x=AML, subset= orig.ident != 'MUTZ3.fresh1')
AML<- subset (x=AML, subset= orig.ident != 'MUTZ3.fresh2')
AML<- subset (x=AML, subset= orig.ident != 'MUTZ3.frozen')
AML<- subset (x=AML, subset= orig.ident != 'OCI.AML3')

AML<- ScaleData(AML)
AML<-FindVariableFeatures(AML)
AML <- RunPCA(AML, assay = "RNA", verbose = FALSE)

#subsetting seurat object on one celltype and based on condition (healthy/disease)----------
#this cell type will be used to identify its target genes and focus on them in the analysis
#so any detected external signals would be targetting this specific cell type
all.samples<-AML
healthy_samples<- subset(x=AML, subset= PredictionRefined == 'normal')
AML_samples<- subset(x=AML, subset= PredictionRefined == 'malignant')

#HSC
HSC.AML<-subset(x= AML_samples, subset= CellType == c('HSC', 'HSC-like')
HSC.AML<-ScaleData(HSC.AML)
HSC.AML<-FindVariableFeatures(HSC.AML)
HSC.AML <- RunPCA(HSC.AML, assay = "RNA", verbose = FALSE)
HSC.AML<- RunUMAP(HSC.AML, dims = 1:20)
DimPlot(HSC.AML, reduction = "umap",  group.by = 'CellType', label= TRUE, label.size=3) +
  ggtitle('HSC.AML')

#Progenitor cells
Prog.AML<-subset(x=AML_samples,subset= CellType == 'Prog-like')

#GMP cells
GMP.AML<- subset(x=AML_samples,subset= CellType == 'GMP-like')

#ProMonocyte cells
ProMono.AML<- subset(x=AML_samples,subset= CellType == 'ProMono-like')

#Monocyte
Mono.AML<- subset(x=AML_samples,subset= CellType == 'Mono-like')

#cDC
cDC.AML<-  subset(x=AML_samples,subset= CellType == 'cDC-like')

#perform exFINDER ---------------

# discovering marker genes------
cDC.AML.markers <- FindAllMarkers(cDC.AML)
cDC.AML.markers<- cDC.AML.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1 & p_val_adj< 0.01)
cDC.AML.markers<- cDC.AML.markers%>% arrange(avg_log2FC,)

#load exFINDER Database-----
load("C:\\Users\\Dell\\Downloads\\bachelor\\exFINDER_method\\data\\LR_layer1_human.rda")
load("C:\\Users\\Dell\\Downloads\\bachelor\\exFINDER_method\\data\\RTF_layer2_human.rda")
load("C:\\Users\\Dell\\Downloads\\bachelor\\exFINDER_method\\data\\TFT_layer3_human.rda")

exFINDER.H <- list()
exFINDER.H[[1]] <- LR_layer1_human
exFINDER.H[[2]] <- RTF_layer2_human
exFINDER.H[[3]] <- TFT_layer3_human

#load CelChat db-----
load("C:\\Users\\Dell\\Downloads\\bachelor\\exFINDER_method\\data\\interaction_input_CellChatDB_human.rda")

#create exFINDER inputs-----

all.samples@assays$RNA$data-> expMydata #input normalized data
row.names(head(expMydata))  #genes should be in rows with rownames and cells in columns with colnames
colnames(expMydata)
meta.AML<- all.samples@meta.data

meta.AML$cell_barcode <- rownames(meta.AML)

meta.AML$Type <- as.factor(meta.AML$CellType)
meta.AML$orig.ident <- as.character(meta.AML$orig.ident)

#check distribution of gene expression------

percentile_data <- get_percentile(Exp.Data = expMydata,
                                  Meta.Data = meta.AML,
                                  percentile = c(.5, .75, .90))

ggplot(data=percentile_data, aes(x=Type, y=Ave.Exp., color=Prob.)) +
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) +
  ggtitle('percentile.AML')

#set up target genes and infer ltGRN------
target.gene <- cDC.AML.markers$gene[1:10]
ltGRN <- get_ltGRN(Target = target.gene, DB = exFINDER.H)  # run everything each time you change no. of target genes


# infer potential external signals from ltGRN----

Pex <- get_potentialex(Graph = ltGRN,
                       Exp.Data = expMydata,
                       Meta.Data = meta.AML,
                       cutoff = -0.1,
                       AG.L = NULL)
head(Pex)

#infer exSigNet----
exSigNet <- get_exSigNet(Graph = ltGRN,
                         Ligands = Pex,
                         Exp.Data = expMydata,
                         Meta.Data = meta.AML,
                         par.cutoff = c(0.5,0.9), 
                         AG.R = c("cDC-like"),
                         AG.TF = c("cDC-like"))


#However, some receptors may interact with both external (come from the external environment) and internal signals (expressed by measured cells), 
#which are challenging to evaluate the effect of external signals.
#We can run the following function to select the receptors that only interact with external signals.

# based on this "exSigNet.0", we identify the receptors that only activated by these external signals 

exSigNet.0 <- filterLR_exSigNet(Graph = exSigNet,
                                Exp.Data = expMydata,
                                DB = exFINDER.H)

#check correlation of number of target genes and transcription factors------

targets.10<- length(exSigNet.0$Graph.Node$Role[exSigNet.0$Graph.Node$Role== 'Target'])
targets.25<- length(exSigNet.0$Graph.Node$Role[exSigNet.0$Graph.Node$Role== 'Target'])
targets.50<- length(exSigNet.0$Graph.Node$Role[exSigNet.0$Graph.Node$Role== 'Target'])
targets.100 <- length(exSigNet.0$Graph.Node$Role[exSigNet.0$Graph.Node$Role== 'Target'])
targets<- as.vector(c(targets.10,targets.25,targets.50,targets.100))


TF.10 <-length(exSigNet.0$Graph.Node$Role[exSigNet.0$Graph.Node$Role== 'TF'])
TF.25 <- length(exSigNet.0$Graph.Node$Role[exSigNet.0$Graph.Node$Role== 'TF'])
TF.50 <- length(exSigNet.0$Graph.Node$Role[exSigNet.0$Graph.Node$Role== 'TF'])
TF.100<- length(exSigNet.0$Graph.Node$Role[exSigNet.0$Graph.Node$Role== 'TF'])
TFs<-as.vector(c(TF.10,TF.25,TF.50,TF.100))

number<-as.vector(c(10,25,50,100))

TF_target<-as.data.frame(targets,TFs)
TF_target$TF<-TFs
rownames(TF_target)<-number

plot(x=number,y= TFs, main = "#targets x #TFs", 
     xlab = "Targets", ylab = "TFs", pch = 19, col = "blue")
lines(x=number, y= TFs, type = "o", col = "black")

# check the inferred external signals and their expressions via heatmap

exS2<- exSigNet.0$Graph.Node$Gene
exS2_exp <- get_AveExp(Gene = exS2,
                      Exp.Data = expMydata,
                      Meta.Data = meta.AML)
exS2_exp.HP <- exS2_exp
exS2_exp.HP <- exS2_exp.HP[, -1]
rownamesHP<- exS2_exp$Gene
row.names(exS2_exp.HP) <- rownamesHP

#if duplicated values are present:
exS2_exp<- exS2_exp[-c(22,37),] # removing duplicate rows
exS2_exp.HP<- exS2_exp.HP [-c(22,37),] #removing duplicated rows
rownamesHP<- exS2_exp$Gene
row.names(exS2_exp.HP) <- rownamesHP
#---------
exS2_exp.HP.mat = as.matrix(exS2_exp.HP) #input has to be matrix
annotation_row = data.frame(Role = exSigNet.0$Graph.Node$Role) 

#if duplicated values are present:
annotation_row <- annotation_row [-c(22,37),]
annotation_row<- as.data.frame(annotation_row)
colnames(annotation_row)<-'Role'

pdf(file ="heatmap.pdf", width =12, height=20)

pheatmap(mat = exS2_exp.HP.mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         cellwidth = 15,
         cellheight = 10,
         gaps_col = NULL,
         gaps_row = c(14,21,40),
         fontsize_row = 10,
         fontsize_col = 10,
         angle_col = NULL,
         legend_breaks = NULL,
         legend_labels = NULL,
         annotation_row = annotation_row) +
         ggtitle('cDC.AML.100')

dev.off()


#Since the ligands come from the external environment, and for the need of quantitative analysis, we set their expression levels equal to 1.

# calculate the expression levels
exSigNet.0 <- get_NodeValue(Graph = exSigNet.0,
                            Exp.Data = expMydata,
                            Meta.Data = meta.AML,
                            AG.R = c('cDC-like'),
                            AG.TF = c('cDC-like'),
                            AG.T = c('cDC-like'),
                            # if there is no expresison data for ligands, exFINDER will set them to 1.
                            Exp.ExData = NULL,
                            Meta.ExData = NULL,
                            AG.ExData = NULL)

# predict signaling strengths
exSigNet.0 <- get_EdgeWeight(Graph = exSigNet.0, 
                             Kh = 2)

#Visualize-----

#hierarchy plot 
exSigNet_plot <- get_readyforplot(Graph = exSigNet.0)

pdf(file ="circleplot.pdf", width =12, height=20)
ggraph(exSigNet_plot[[1]], layout = 'dendrogram', circular = TRUE) + 
  geom_node_point(aes(filter = leaf, x=x, y=y, colour=Role, size=Expression)) +
  theme_void() +
  geom_conn_bundle(data = get_con(from = exSigNet_plot[[2]]$from, 
                                  to = exSigNet_plot[[2]]$to,
                                  Strength = exSigNet_plot[[2]]$weight), width = 1.0, aes(colour=Strength),
                   arrow = arrow(length = unit(3, 'mm')), start_cap = circle(3, 'mm'), end_cap = circle(3, 'mm')) +
  scale_edge_color_continuous(low="#80b1d3", high="orange")+
  geom_node_text(aes(label = node)) +
  scale_size_continuous( range = c(0.1,10) ) +
  ggtitle('Circle Plot cDC.AML.100')
dev.off()
#hive plot
Graph <- exSigNet.0

Role <- factor(x = exSigNet.0$Graph.Node$Role, levels = c("Ligand", "Receptor", "TF", "Target"))
edges <- data.frame("from" = exSigNet.0$Graph.Edge$from, "to" = exSigNet.0$Graph.Edge$to)

edges <- edges %>%
  mutate(weight = exSigNet.0$Graph.Edge$Weight)

nodes <- data.frame(name = exSigNet.0$Graph.Node$Gene, role = Role, value = exSigNet.0$Graph.Node$Size)
g2 <- graph_from_data_frame(vertices = nodes, d = edges)

pdf(file ="hiveplot.pdf", width =12, height=20)
ggraph(g2, 'hive', axis = nodes$role) + 
  geom_edge_hive(aes(colour = weight), width = 0.7, arrow = arrow(length = unit(3, 'mm')), start_cap = circle(3, 'mm'), end_cap = circle(3, 'mm')) + 
  # geom_axis_hive(aes(colour = role), size = 0, label = FALSE) + 
  geom_node_point(aes(size = value, color = role),  alpha = 1) +
  geom_node_text(aes(label = name), angle = 45,  size = 3, repel = TRUE) +
  coord_fixed() +
  scale_edge_color_continuous(low="#80b1d3", high="#FF7F24") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.ticks =element_blank(),
    axis.text =element_blank(),
    axis.title = element_blank()
  ) +
  ggtitle('cDC.AML.100')
dev.off()

#Predict critical external signals and targets based on total signal flows-----
df <- get_LRMaxFlow(exSigNet.0)
df.TotalFlow <- compare_MaxFlow(df = df)
df.TotalFlow 

df.TotalFlow %>%
  mutate(Gene = fct_reorder(Gene, TotalFlow* (Role=='Ligand'))) %>%
  mutate(Gene = fct_reorder(Gene, TotalFlow* (Role=='Target'))) %>%
  ggplot(aes(x=Gene, y=TotalFlow,  fill=Role)) +
  geom_bar(stat="identity", alpha=1, width=.6) +
  xlab("Gene") +
  ylab("Total signal flow") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))+
  ggtitle('Total Signal Flow HSC.AML')

#Clustering signaling networks between ligand-target pairs------

sm<-get_SimilarityMatrix.0(Graph = exSigNet.0)
desc_stats = data.frame( Min=apply(sm, 2, min),#minimum
                         Med=apply(sm, 2, median),#median
                         Mean=apply(sm, 2, mean),#mean
                         SD=apply(sm, 2, sd),#Standard deviation
                         Max=apply(sm, 2, max)#maximum
)
desc_stats = round(desc_stats, 1)

library(factoextra)
library(cluster)
sm = scale(sm)
sm = sm[, colSums(is.na(sm)) != nrow(sm)]

res = get_clust_tendency(sm, 10, graph = TRUE)
km.res = kmeans(sm, 4, nstart = 20)

ligand_target_pairs<- rownames(sm)
fviz_cluster(km.res, data = sm, geom = "point", 
             label = TRUE,
             pointsize = 3.5)+
  theme_minimal()+
  geom_text(aes(label = ligand_target_pairs, size = 0.2, vjust = -1))

 
 