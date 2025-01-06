# Project: detecting external signals that can drive disease and fibrosis in AML Laila Mohamed, 2024

# load required libraries

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

# Download data-------

proj_dir <- getwd()
AML <- readRDS(file.path(proj_dir, "data_whole_processed/Seurat_AML_sparse.rds"))

# checking presence of all samples and celltypes-------
unique(AML@meta.data$CellType)
unique(AML@meta.data$orig.ident)

# (t(table(AML@meta.data$orig.ident, AML@meta.data$CellType)) == 0) %>% lattice::levelplot()

((table(AML@meta.data$PredictionRefined, AML@meta.data$CellType))) %>% lattice::levelplot()
AML@meta.data %>% as_tibble()
AML@meta.data$orig.ident %>% table


# quality control and umaps------
AML <- subset(x = AML, subset = orig.ident != "MUTZ3.fresh1")
AML <- subset(x = AML, subset = orig.ident != "MUTZ3.fresh2")
AML <- subset(x = AML, subset = orig.ident != "MUTZ3.frozen")
AML <- subset(x = AML, subset = orig.ident != "OCI.AML3")

all.samples <- AML
rm(AML) # to save space

healthy_samples <- subset(x = all.samples, 
                          subset = orig.ident %in% 
                            c("BM1", "BM2", "BM3", "BM4", "BM5.34p",
                              "BM5.34p38n"))

AML_samples <- subset(
  x = all.samples,
  subset = orig.ident %in% 
    c("AML1012.D0", "AML210A.D0", "AML419A.D0", "AML916.D0",
      "AML921A.D0", "AML314.D0", "AML314.D31", "AML371.D0", "AML371.D34", 
      "AML475.D0", "AML475.D29", "AML722B.D0",
      "AML722B.D49", "AML870.D0", " AML870.D14", " AML997.D0", "AML997.D35",
      "AML329.D0", "AML329.D20", "AML329.D37",
      "AML420B.D0", "AML420B.D14 ", "AML420B.D35", "AML556.D0", "AML556.D15",
      "AML556.D31", "AML328.D0", "AML328.D29 ",
      "AML328.D113", "AML328.D171", "AML707B.D0", "AML707B.D18", "AML707B.D41",
      "AML707B.D97", "AML707B.D113"))

malignant_cells <- subset(x = AML_samples, subset = PredictionRefined == "malignant")
healthycells_patient <- subset(x = AML_samples, subset = PredictionRefined == "normal")

malignant_cells <- ScaleData(malignant_cells)
malignant_cells <- FindVariableFeatures(malignant_cells)
malignant_cells <- RunPCA(malignant_cells, assay = "RNA", verbose = FALSE)
malignant_cells <- RunUMAP(malignant_cells, dims = 1:20)
getOption("formatR.width", getOption("width"))

DimPlot(malignant_cells, reduction = "umap", group.by = "CellType", label = TRUE, label.size = 3) + ggtitle("Cell Types Malignant")

donor_malignant <- subset(x = all.samples, subset = orig.ident %in% c("BM1", "BM2", "BM3", "BM4", "BM5.34p",
                                                                      "BM5.34p38n") | PredictionRefined == "malignant")

AML_samples <- ScaleData(AML_samples)
AML_samples <- FindVariableFeatures(AML_samples)  #I want to detect the marker genes that differentiate between malignant and nonmalignant cells in AML sample
AML_samples <- RunPCA(AML_samples, assay = "RNA", verbose = FALSE)
AML_samples <- RunUMAP(AML_samples, dims = 1:20)

DimPlot(AML_samples, reduction = "umap", group.by = "PredictionRefined",
        label = TRUE, label.size = 3) + ggtitle("AML_samples state")

healthy_samples <- ScaleData(healthy_samples)
healthy_samples <- FindVariableFeatures(healthy_samples)  #I want to detect the marker genes that differentiate between malignant and nonmalignant cells in AML sample and cell types
healthy_samples <- RunPCA(healthy_samples, assay = "RNA", verbose = FALSE)
healthy_samples <- RunUMAP(healthy_samples, dims = 1:20)

## remove unused factors
healthy_samples@meta.data = healthy_samples@meta.data %>% 
  mutate(CellType = fct_drop(CellType)) %>% 
  mutate(PredictionRefined = fct_drop(PredictionRefined))

DimPlot(healthy_samples, reduction = "umap", group.by = "PredictionRefined", 
        label = TRUE, label.size = 3) +
  ggtitle("healthy_samples state")

healthy_samples@meta.data$CellType %>% table
DimPlot(healthy_samples, reduction = "umap", group.by = "CellType", 
        label = TRUE, label.size = 3) +
  ggtitle("healthy_samples state")

donor_malignant <- ScaleData(donor_malignant)
donor_malignant <- FindVariableFeatures(donor_malignant)  #I want to detect the marker genes that differentiate between conditions
donor_malignant <- RunPCA(donor_malignant, assay = "RNA", verbose = FALSE)
donor_malignant <- RunUMAP(donor_malignant, dims = 1:20)

DimPlot(donor_malignant, reduction = "umap", group.by = "PredictionRefined",
        label = TRUE, label.size = 3) +
  ggtitle("healthy_malignant state")

DimPlot(donor_malignant, reduction = "umap", group.by = "CellType",
        label = TRUE, label.size = 3) +
  ggtitle("healthy_malignant state")

healthycells_patient <- ScaleData(healthycells_patient)
healthycells_patient <- FindVariableFeatures(healthycells_patient)  #I want to detect the marker genes that differentiate between malignant and nonmalignant cells in AML sample
healthycells_patient <- RunPCA(healthycells_patient, assay = "RNA", verbose = FALSE)
healthycells_patient <- RunUMAP(healthycells_patient, dims = 1:20)




# subsetting seurat object on one celltype and based on condition (healthy/disease)----------this cell
# type will be used to identify its target genes and focus on them in the analysis so any detected external
# signals would be targetting this specific cell type rhis step seems to be unnecessary

# HSC
HSC.AML <- subset(x = AML_samples, subset = CellType == c("HSC", "HSC-like"))
HSC.AML <- ScaleData(HSC.AML)
HSC.AML <- FindVariableFeatures(HSC.AML)
HSC.AML <- RunPCA(HSC.AML, assay = "RNA", verbose = FALSE)
HSC.AML <- RunUMAP(HSC.AML, dims = 1:20)
DimPlot(HSC.AML, reduction = "umap", group.by = "CellType",
        label = TRUE, label.size = 3) + ggtitle("HSC.AML")

# Progenitor cells
Prog.Patient <- subset(x = AML_samples, subset = CellType == c("Prog", "Prog-like"))
Prog.malignant <- subset(x = malignant_cells, subset = CellType == "Prog-like")  #GMP cells
GMP.PatientL <- subset(x = AML_samples, subset = CellType == c("GMP", "GMP-like"))
GMP.AML <- ScaleData(GMP.AML)

# ProMonocyte cells
ProMono.AML <- subset(x = AML_samples, subset = CellType == "ProMono-like")
ProMono.AML <- ScaleData(ProMono.AML)

# Monocyte
Mono.AML <- subset(x = AML_samples, subset = CellType == "Mono-like")
Mono.AML <- ScaleData(Mono.AML)

# cDC
cDC.AML <- subset(x = AML_samples, subset = CellType == "cDC-like")
cDC.AML <- ScaleData(cDC.AML)
cDC.AML <- FindVariableFeatures(cDC.AML)
cDC.AML <- RunPCA(cDC.AML, assay = "RNA", verbose = FALSE)
cDC.AML <- RunUMAP(cDC.AML, dims = 1:20)
DimPlot(HSC.AML, reduction = "umap", group.by = "CellType", label = TRUE, label.size = 3) + ggtitle("HSC.AML")

# discovering marker genes------
Idents(healthy_samples) = "CellType"
healthy.markers <- FindAllMarkers(healthy_samples)
healthy.markers <- healthy.markers %>%
  dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.01)

## Option1: consider the top markers by avg_log2FC
markers_by_log2fc <- healthy.markers %>%
  as_tibble() %>% 
  # group_by(cluster) %>% 
  arrange(desc(avg_log2FC))

VlnPlot(healthy_samples, markers_by_log2fc$gene[1:6]) + 
  ggtitle("Top X Marker Genes by log2FC")

## Option2: consider the top markers by p_val_adj
markers_by_p_val = markers_by_log2fc %>% 
  arrange(p_val_adj)

VlnPlot(healthy_samples, markers_by_p_val$gene[1:12]) &
  patchwork::plot_annotation(title = "Top X Marker Genes by p_val")
ggsave("markers_violin_HH_top_p_val.pdf", height = 20, width = 30,
       path = file.path(proj_dir, "results/marker_genes"))

FeaturePlot(healthy_samples, features = markers_by_p_val$gene[1:12]) + 
  ggtitle("Top X Marker Genes by p_val")
ggsave("markers_featureplot_HH_top_p_val.pdf", height = 20, width = 30,
       path = file.path(proj_dir, "results/marker_genes"))

## Plot top1 marker gene for each cluster
library(magrittr)
top1genes_p_val = markers_by_p_val %>% group_by(cluster) %>% arrange(p_val_adj) %>% 
  slice_head(n = 1) %>% ungroup()
FeaturePlot(healthy_samples, features = top1genes_p_val$gene) +
  plot_annotation(tag_levels = list(top1genes_p_val$cluster), 
                  title = "Top1 Marker Genes (by p_val) for each CellType")
ggsave("markers_featureplot_HH_top1_p_val.pdf", height = 20, width = 30,
       path = file.path(proj_dir, "results/marker_genes"))

## Plot marker genes heatmaps
top_markers_by_p_val = markers_by_p_val %>%
  group_by(cluster) %>%
  slice_head(n = 10) %>%
  ungroup()
top_markers_by_log2fc = markers_by_log2fc %>%
  group_by(cluster) %>%
  slice_head(n = 10) %>%
  ungroup()

DoHeatmap(healthy_samples, features = top_markers_by_p_val$gene) + NoLegend()
ggsave("markers_heatmap_HH_top_p_val.pdf", height = 20, width = 30,
       path = file.path(proj_dir, "results/marker_genes"))

DoHeatmap(healthy_samples, features = top_markers_by_log2fc$gene) + NoLegend()
ggsave("markers_heatmap_HH_top_log2fc.pdf", height = 20, width = 30,
       path = file.path(proj_dir, "results/marker_genes"))



malignant.markers <- c("S100A9", "CD74", "CTSS", "EMB", "PSAP", "PCNP", "CFD", "APLP2", "EEF1A1", "RPS3A", "VIM",
                       "HSPA5", "NPM1", "NEAT1", "HMGB2", "HLA-DRB1", "CD74", "TXNIP", "EEF2", "ANXA1")
#'RPL13','RPL32','RPL41','RPL7','RPL9','RPS13','RPS19','RPS21','RPS29','RPS4X','RPS14','RPS7','RPSA'

target.gene.healthy <- healthy.markers$gene[84:95]
target.gene <- c("CTC", "HBBP1", "HBET", "KLRC3", "PRF1", "CD3G", "FGR", "NRIP1", "TYROBP", "MME", "BCL11B",
                 "RASSF4", "SAP30", "FGL2")
DoHeatmap(healthy_samples, features = target.gene.healthy, group.by = "CellType") + ggtitle("Healthy Marker Genes log2FC>1")
FeaturePlot(healthy_samples, features = target.gene.healthy)

# perform exFINDER ---------------

# load exFINDER Database-----
load("C:\\Users\\Dell\\Downloads\\bachelor\\exFINDER_method\\data\\LR_layer1_human.rda")
load("C:\\Users\\Dell\\Downloads\\bachelor\\exFINDER_method\\data\\RTF_layer2_human.rda")
load("C:\\Users\\Dell\\Downloads\\bachelor\\exFINDER_method\\data\\TFT_layer3_human.rda")

exFINDER.H <- list()
exFINDER.H[[1]] <- LR_layer1_human
exFINDER.H[[2]] <- RTF_layer2_human
exFINDER.H[[3]] <- TFT_layer3_human

# load CelChat db-----
load("C:\\Users\\Dell\\Downloads\\bachelor\\exFINDER_method\\data\\interaction_input_CellChatDB_human.rda")

# create exFINDER inputs-----

all.samples@assays$RNA$data -> expMydata  #input normalized data
row.names(head(expMydata))  #genes should be in rows with rownames and cells in columns with colnames
colnames(expMydata)
meta.AML <- all.samples@meta.data

meta.AML$cell_barcode <- rownames(meta.AML)

meta.AML$Type <- as.factor(meta.AML$CellType)
meta.AML$orig.ident <- as.character(meta.AML$orig.ident)

# check distribution of gene expression------

percentile_data <- get_percentile(Exp.Data = expMydata, Meta.Data = meta.AML, percentile = c(0.5, 0.75, 0.9))

ggplot(data = percentile_data, aes(x = Type, y = Ave.Exp., color = Prob.)) + geom_point(size = 3) + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) + ggtitle("percentile.AML")

# set up target genes and infer ltGRN------
target.gene <- healthy_malignant.markers$gene[1:100]
ltGRN <- get_ltGRN(Target = target.gene, DB = exFINDER.H)  # run everything each time you change no. of target genes


# infer potential external signals from ltGRN----

Pex <- get_potentialex(Graph = ltGRN, Exp.Data = expMydata, Meta.Data = meta.AML, cutoff = -0.1, AG.L = NULL)
head(Pex)

# infer exSigNet----
exSigNet <- get_exSigNet(Graph = ltGRN, Ligands = Pex, Exp.Data = expMydata, Meta.Data = meta.AML, par.cutoff = c(0.5,
                                                                                                                  0.9), AG.R = c("cDC-like"), AG.TF = c("cDC-like"))


# However, some receptors may interact with both external (come from the external environment) and internal
# signals (expressed by measured cells), which are challenging to evaluate the effect of external signals.
# We can run the following function to select the receptors that only interact with external signals.

# based on this 'exSigNet.0', we identify the receptors that only activated by these external signals

exSigNet.0 <- filterLR_exSigNet(Graph = exSigNet, Exp.Data = expMydata, DB = exFINDER.H)

# check correlation of number of target genes and transcription factors------

targets.10 <- length(exSigNet.0$Graph.Node$Role[exSigNet.0$Graph.Node$Role == "Target"])
targets.25 <- length(exSigNet.0$Graph.Node$Role[exSigNet.0$Graph.Node$Role == "Target"])
targets.50 <- length(exSigNet.0$Graph.Node$Role[exSigNet.0$Graph.Node$Role == "Target"])
targets.100 <- length(exSigNet.0$Graph.Node$Role[exSigNet.0$Graph.Node$Role == "Target"])
targets <- as.vector(c(targets.10, targets.25, targets.50, targets.100))


TF.10 <- length(exSigNet.0$Graph.Node$Role[exSigNet.0$Graph.Node$Role == "TF"])
TF.25 <- length(exSigNet.0$Graph.Node$Role[exSigNet.0$Graph.Node$Role == "TF"])
TF.50 <- length(exSigNet.0$Graph.Node$Role[exSigNet.0$Graph.Node$Role == "TF"])
TF.100 <- length(exSigNet.0$Graph.Node$Role[exSigNet.0$Graph.Node$Role == "TF"])
TFs <- as.vector(c(TF.10, TF.25, TF.50, TF.100))

number <- as.vector(c(10, 25, 50, 100))

TF_target <- as.data.frame(targets, TFs)
TF_target$TF <- TFs
rownames(TF_target) <- number

plot(x = number, y = TFs, main = "#targets x #TFs", xlab = "Targets", ylab = "TFs", pch = 19, col = "blue")
lines(x = number, y = TFs, type = "o", col = "black")

# check the inferred external signals and their expressions via heatmap

exS2 <- exSigNet.0$Graph.Node$Gene
exS2_exp <- get_AveExp(Gene = exS2, Exp.Data = expMydata, Meta.Data = meta.AML)
exS2_exp.HP <- exS2_exp
exS2_exp.HP <- exS2_exp.HP[, -1]
rownamesHP <- exS2_exp$Gene
row.names(exS2_exp.HP) <- rownamesHP

# if duplicated values are present:
exS2_exp <- exS2_exp[-c(22, 37), ]  # removing duplicate rows
exS2_exp.HP <- exS2_exp.HP[-c(22, 37), ]  #removing duplicated rows
rownamesHP <- exS2_exp$Gene
row.names(exS2_exp.HP) <- rownamesHP
#---------
exS2_exp.HP.mat <- as.matrix(exS2_exp.HP)  #input has to be matrix
annotation_row <- data.frame(Role = exSigNet.0$Graph.Node$Role)

# if duplicated values are present:
annotation_row <- annotation_row[-c(22, 37), ]
annotation_row <- as.data.frame(annotation_row)
colnames(annotation_row) <- "Role"

pdf(file = "heatmap.pdf", width = 12, height = 20)

pheatmap(mat = exS2_exp.HP.mat, cluster_rows = FALSE, cluster_cols = FALSE, cellwidth = 15, cellheight = 10,
         gaps_col = NULL, gaps_row = c(14, 21, 40), fontsize_row = 10, fontsize_col = 10, angle_col = NULL, legend_breaks = NULL,
         legend_labels = NULL, annotation_row = annotation_row) + ggtitle("cDC.AML.100")

dev.off()


# Since the ligands come from the external environment, and for the need of quantitative analysis, we set
# their expression levels equal to 1.

# calculate the expression levels
exSigNet.0 <- get_NodeValue(Graph = exSigNet.0, Exp.Data = expMydata, Meta.Data = meta.AML, AG.R = c("cDC-like"),
                            AG.TF = c("cDC-like"), AG.T = c("cDC-like"), # if there is no expresison data for ligands, exFINDER will set them to 1.  
                            Exp.ExData = NULL,  
                            Meta.ExData = NULL,  
                            AG.ExData = NULL))
                
                            
                            # predict signaling strengths
                            exSigNet.0 <- get_EdgeWeight(Graph = exSigNet.0, Kh = 2)
                            
                            # Visualize-----
                            
                            # hierarchy plot
                            
                            # exSigNet_plot <- get_readyforplot(Graph = exSigNet.0)
                            
                            pdf(file = "circleplot.pdf", width = 12, height = 20)
                            ggraph(exSigNet_plot[[1]], layout = "dendrogram", circular = TRUE) + geom_node_point(aes(filter = leaf, x = x,
                                                                                                                     y = y, colour = Role, size = Expression)) + theme_void() + geom_conn_bundle(data = get_con(from = exSigNet_plot[[2]]$from,
                                                                                                                                                                                                                to = exSigNet_plot[[2]]$to, Strength = exSigNet_plot[[2]]$weight), width = 1, aes(colour = Strength), arrow = arrow(length = unit(3,
                                                                                                                                                                                                                                                                                                                                                  "mm")), start_cap = circle(3, "mm"), end_cap = circle(3, "mm")) + scale_edge_color_continuous(low = "#80b1d3",
                                                                                                                                                                                                                                                                                                                                                                                                                                                high = "orange") + geom_node_text(aes(label = node)) + scale_size_continuous(range = c(0.1, 10)) + ggtitle("Circle Plot cDC.AML.100")
                            dev.off()
                            # hive plot
                            Graph <- exSigNet.0
                            
                            Role <- factor(x = exSigNet.0$Graph.Node$Role, levels = c("Ligand", "Receptor", "TF", "Target"))
                            edges <- data.frame(from = exSigNet.0$Graph.Edge$from, to = exSigNet.0$Graph.Edge$to)
                            
                            edges <- edges %>%
                              mutate(weight = exSigNet.0$Graph.Edge$Weight)
                            
                            nodes <- data.frame(name = exSigNet.0$Graph.Node$Gene, role = Role, value = exSigNet.0$Graph.Node$Size)
                            g2 <- graph_from_data_frame(vertices = nodes, d = edges)
                            
                            pdf(file = "hiveplot.pdf", width = 12, height = 20)
                            ggraph(g2, "hive", axis = nodes$role) + geom_edge_hive(aes(colour = weight), width = 0.7, arrow = arrow(length = unit(3,
                                                                                                                                                  "mm")), start_cap = circle(3, "mm"), end_cap = circle(3, "mm")) + # geom_axis_hive(aes(colour = role), size = 0, label = FALSE) +  "mm")),
                              "mm")), start_cap = circle(3, "mm"), end_cap = circle(3, "mm")) + # geom_axis_hive(aes(colour = role), size = 0, label = FALSE) +  start_cap
  "mm")), start_cap = circle(3, "mm"), end_cap = circle(3, "mm")) + # geom_axis_hive(aes(colour = role), size = 0, label = FALSE) +  =
  "mm")), start_cap = circle(3, "mm"), end_cap = circle(3, "mm")) + # geom_axis_hive(aes(colour = role), size = 0, label = FALSE) +  circle(3,
  "mm")), start_cap = circle(3, "mm"), end_cap = circle(3, "mm")) + # geom_axis_hive(aes(colour = role), size = 0, label = FALSE) +  "mm"),
  "mm")), start_cap = circle(3, "mm"), end_cap = circle(3, "mm")) + # geom_axis_hive(aes(colour = role), size = 0, label = FALSE) +  end_cap
  "mm")), start_cap = circle(3, "mm"), end_cap = circle(3, "mm")) + # geom_axis_hive(aes(colour = role), size = 0, label = FALSE) +  =
  "mm")), start_cap = circle(3, "mm"), end_cap = circle(3, "mm")) + # geom_axis_hive(aes(colour = role), size = 0, label = FALSE) +  circle(3,
  "mm")), start_cap = circle(3, "mm"), end_cap = circle(3, "mm")) + # geom_axis_hive(aes(colour = role), size = 0, label = FALSE) +  "mm"))
  "mm")), start_cap = circle(3, "mm"), end_cap = circle(3, "mm")) + # geom_axis_hive(aes(colour = role), size = 0, label = FALSE) +  +
  "mm")), start_cap = circle(3, "mm"), end_cap = circle(3, "mm")) + # geom_axis_hive(aes(colour = role), size = 0, label = FALSE) +  #
  "mm")), start_cap = circle(3, "mm"), end_cap = circle(3, "mm")) + # geom_axis_hive(aes(colour = role), size = 0, label = FALSE) +  geom_axis_hive(aes(colour
  "mm")), start_cap = circle(3, "mm"), end_cap = circle(3, "mm")) + # geom_axis_hive(aes(colour = role), size = 0, label = FALSE) +  =
  "mm")), start_cap = circle(3, "mm"), end_cap = circle(3, "mm")) + # geom_axis_hive(aes(colour = role), size = 0, label = FALSE) +  role),
  "mm")), start_cap = circle(3, "mm"), end_cap = circle(3, "mm")) + # geom_axis_hive(aes(colour = role), size = 0, label = FALSE) +  size
  "mm")), start_cap = circle(3, "mm"), end_cap = circle(3, "mm")) + # geom_axis_hive(aes(colour = role), size = 0, label = FALSE) +  =
  "mm")), start_cap = circle(3, "mm"), end_cap = circle(3, "mm")) + # geom_axis_hive(aes(colour = role), size = 0, label = FALSE) +  0,
  "mm")), start_cap = circle(3, "mm"), end_cap = circle(3, "mm")) + # geom_axis_hive(aes(colour = role), size = 0, label = FALSE) +  label
  "mm")), start_cap = circle(3, "mm"), end_cap = circle(3, "mm")) + # geom_axis_hive(aes(colour = role), size = 0, label = FALSE) +  =
  "mm")), start_cap = circle(3, "mm"), end_cap = circle(3, "mm")) + # geom_axis_hive(aes(colour = role), size = 0, label = FALSE) +  FALSE)
  "mm")), start_cap = circle(3, "mm"), end_cap = circle(3, "mm")) + # geom_axis_hive(aes(colour = role), size = 0, label = FALSE) +  +
  geom_node_point(aes(size = value, color = role), alpha = 1) + geom_node_text(aes(label = name), angle = 45, size = 3,
                                                                               repel = TRUE) + coord_fixed() + scale_edge_color_continuous(low = "#80b1d3", high = "#FF7F24") + theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
        axis.title = element_blank()) + ggtitle("cDC.AML.100")
dev.off()

# Predict critical external signals and targets based on total signal flows-----
df <- get_LRMaxFlow(exSigNet.0)
df.TotalFlow <- compare_MaxFlow(df = df)
df.TotalFlow

df.TotalFlow %>%
  mutate(Gene = fct_reorder(Gene, TotalFlow * (Role == "Ligand"))) %>%
  mutate(Gene = fct_reorder(Gene, TotalFlow * (Role == "Target"))) %>%
  ggplot(aes(x = Gene, y = TotalFlow, fill = Role)) + geom_bar(stat = "identity", alpha = 1, width = 0.6) +
  xlab("Gene") + ylab("Total signal flow") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 0.5,
                                                                                           vjust = 0.5)) + ggtitle("Total Signal Flow HSC.AML")

# Clustering signaling networks between ligand-target pairs------ 

sm<-get_SimilarityMatrix.0(Graph = exSigNet.0)
desc_stats = data.frame( Min=apply(sm, 2, min), # minimum
                         Med=apply(sm, 2, median), # median
                         Mean=apply(sm, 2, mean), # mean
                         SD=apply(sm, 2, sd), # Standard deviation
                         Max=apply(sm, 2, max) # maximum
)

desc_stats <- round(desc_stats, 1)

library(factoextra)
library(cluster)
sm <- scale(sm)
sm <- sm[, colSums(is.na(sm)) != nrow(sm)]

res <- get_clust_tendency(sm, 10, graph = TRUE)
km.res <- kmeans(sm, 4, nstart = 20)

ligand_target_pairs <- rownames(sm)
fviz_cluster(km.res, data = sm, geom = "point", label = TRUE, pointsize = 3.5) + theme_minimal() + geom_text(aes(label = ligand_target_pairs,
                                                                                                                 size = 0.2, vjust = -1))
