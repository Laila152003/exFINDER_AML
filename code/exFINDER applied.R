# applying exFINDER on the mouse dataset from https://doi.org/10.1016/j.stem.2020.11.004
#installing exFINDER----
devtools::install_github("ChanghanGitHub/exFINDER")
 #load libraries-----
library(exFINDER)
library(Seurat)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(dplyr)
library(ComplexHeatmap)
library(ggraph)
library(stringr)
library(igraph)
#read data----
mydata<-readRDS("C:/Users/Dell/Downloads/bachelor/jak2_tpo.rds")
mydata@meta.data$Type<- as.factor(mydata@meta.data$anno)
#checking cell number in each celltype---------
cellcounts<-table(mydata@meta.data$anno)
#adipogenic MSC: 1233 #osteogenic MSC: 1040    #transition MSC: 478
# interferon high MSCs: 189     #nonmyelinating schwann cell progenitors: 170  #myelinating SCPs: 183
#Osteoblastic lineage cells: 184  #endothelial cells: 109

# discovering marker genes------

MF.markers <- FindAllMarkers(mydata)
MF.markers<- MF.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1 & p_val_adj< 0.01)
MF.marker<- MF.markers%>% arrange(avg_log2FC)


#load exFINDER Database-----
load("C:\\Users\\Dell\\Downloads\\bachelor\\exFINDER\\data\\LR_layer1_mouse.rda")
load("C:\\Users\\Dell\\Downloads\\bachelor\\exFINDER\\data\\RTF_layer2_mouse.rda")
load("C:\\Users\\Dell\\Downloads\\bachelor\\exFINDER\\data\\TFT_layer3_mouse.rda")

exFINDER.M <- list()
exFINDER.M[[1]] <- LR_layer1_mouse
exFINDER.M[[2]] <- RTF_layer2_mouse
exFINDER.M[[3]] <- TFT_layer3_mouse

#load CelChat db-----
load("C:\\Users\\Dell\\Downloads\\bachelor\\exFINDER\\data\\interaction_input_CellChatDB_mouse.rda")

#create exFINDER inputs-----
mydata@assays$RNA$data-> expMydata 
row.names(expMydata)  #genes should be in rows with rownames and cells in columns with colnames
colnames(mydata)
meta.mydata<- mydata@meta.data

meta$cell_barcode <- rownames(meta)
colnames(Data) <- t(meta$cell_barcode)

meta$Type <- as.character(meta$Type)
meta$orig.ident <- as.character(meta$orig.ident)


#Infer ltGRN------------

#set up target genes and infer ltGRN------
target.gene <- MF.marker$gene[1:2000]
ltGRN <- get_ltGRN(Target = target.gene, DB = exFINDER.M)


percentile_data <- get_percentile(Exp.Data = expMydata,
                                  Meta.Data = meta.mydata,
                                  percentile = c(.5, .75, .90))

pgraph<- ggplot(data=percentile_data, aes(x=Type, y=Ave.Exp., color=Prob.)) +
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
pgraph

# infer potential external signals from ltGRN----
Pex <- get_potentialex(Graph = ltGRN,
                       Exp.Data = expMydata,
                       Meta.Data = meta.mydata,
                       cutoff = -0.1,
                       AG.L = NULL)
head(Pex)

#infer exSigNet----
exSigNet <- get_exSigNet(Graph = ltGRN,
                         Ligands = Pex,
                         Exp.Data = expMydata,
                         Meta.Data = meta.mydata,
                         par.cutoff = c(1, 1), 
                         AG.R = c("1, adipogenic MSC", "2, osteogenic MSC"),
                         AG.TF = c("1, adipogenic MSC"))

exSigNet.0 <- filterLR_exSigNet(Graph = exSigNet,
                              Exp.Data = expMydata,
                              DB = exFINDER.M)
However, some receptors may interact with both external (come from the external environment) and internal signals (expressed by measured cells), 
which are challenging to evaluate the effect of external signals.
We can run the following function to select the receptors that only interact with external signals.

# based on this "exSigNet.2", we identify the receptors that only activated by these external signals 
exSigNet.2 <- filterLR_exSigNet(Graph = exSigNet.2,
                  Exp.Data = exp.Exdata,
                  DB = exFINDER.Z)

# check the inferred external signals and their expressions via heatmap
exS.2 <- exSigNet.2$Graph.Node$Gene[exSigNet.2$Graph.Node$Role == "Ligand"]

get_heatmap(Gene = exS.2,
            Exp.Data = exp.Exdata,
            Meta.Data = meta.Exdata)

#Since the ligands come from the external environment, and for the need of quantitative analysis, we set their expression levels equal to 1.

# calculate the expression levels
exSigNet.2 <- get_NodeValue(Graph = exSigNet.2,
                          Exp.Data = exp.Exdata,
                          Meta.Data = meta.Exdata,
                          AG.R = c("skeletal", "transition"),
                          AG.TF = c('skeletal'),
                          AG.T = c('skeletal'),
                          # if there is no expresison data for ligands, exFINDER will set them to 1.
                          Exp.ExData = NULL,
                          Meta.ExData = NULL,
                          AG.ExData = NULL)

# predict signaling strengths
exSigNet.2 <- get_EdgeWeight(Graph = exSigNet.2, 
               Kh = 2)

#Visualize-----
exS <- exSigNet.0$Graph.Node$Gene[exSigNet.0$Graph.Node$Role == "Ligand"]


mat.meta.mydata<-as.matrix(meta.mydata)
get_heatmap(Gene = exS,
            Exp.Data = expMydata,
            Meta.Data =meta.mydata)
#calculate the expression levels----
exSigNet <- get_NodeValue(Graph = exSigNet.0,
                          Exp.Data = expMydata,
                          Meta.Data = meta.mydata,
                          AG.R = c("1, adipogenic MSC", "2, osteogenic MSC"),
                          AG.TF = c("1, adipogenic MSC"),
                          AG.T = c("1, adipogenic MSC"),
                          Exp.ExData = NULL,
                          Meta.ExData = NULL,
                          AG.ExData = NULL)
exSigNet <- get_EdgeWeight(Graph = exSigNet, 
                           Kh = 2)
#visualize--------
exS_exp <- get_AveExp(Gene = exSigNet$Graph.Node$Gene,
                      Exp.Data = expMydata,
                      Meta.Data = meta.mydata)

exS_exp.HP <- exS_exp
exS_exp.HP <- exS_exp.HP[, -1]
row.names(exS_exp.HP) <- exS_exp$Gene
exS_exp.HP->exS_exp
rm(exS_exp.HP)

annotation_row <- data.frame(Role = exSigNet$Graph.Node$Role)
row.names(annotation_row) <- row.names(exS_exp)

#heatmap
pdf(file ="file.pdf", width =12, height=20)
pheatmap(exS_exp,
         cluster_row = FALSE, cluster_cols = FALSE, 
         cellwidth = 12, cellheight = 12,
         gaps_row = c(4, 8, 17),
         fontsize_row = 10,
         fontsize_col = 10,
         annotation_row = annotation_row)
dev.off()
#hierarchy plot 
exSigNet_plot <- get_readyforplot(Graph = exSigNet)

ggraph(exSigNet_plot[[1]], layout = 'dendrogram', circular = TRUE) + 
  geom_node_point(aes(filter = leaf, x=x, y=y, colour=Role, size=Expression)) +
  theme_void() +
  geom_conn_bundle(data = get_con(from = exSigNet_plot[[2]]$from, 
                                  to = exSigNet_plot[[2]]$to,
                                  Strength = exSigNet_plot[[2]]$weight), width = 1.0, aes(colour=Strength),
                   arrow = arrow(length = unit(3, 'mm')), start_cap = circle(3, 'mm'), end_cap = circle(3, 'mm')) +
  scale_edge_color_continuous(low="#80b1d3", high="orange")+
  geom_node_text(aes(label = node)) +
  scale_size_continuous( range = c(0.1,10) )

#Predict critical external signals and targets based on total signal flows-----
df <- get_LRMaxFlow(exSigNet)
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
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))

#Clustering signaling networks between ligand-target pairs------

sm<-get_SimilarityMatrix(Graph = exSigNet)
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
res$plot
km.res = kmeans(sm, 4, nstart = 20)
autoplot(km.res, data=sm, label=TRUE, label.size=4, frame=FALSE) + theme_bw()

#infer functionally related signalling pathway-------
exSigNet.Path <- check_Pathway(Graph = exSigNet,
                               LR.DB = CellChatLR.mouse,
                               Par = 2)

exSigNet.Path #NA? :/

#Evaluate GO analysis results------
library(org.Mm.eg.db)
GO.exSigNet <- get_EnrichAnalysis(Graph = exSigNet,
                                  OrgDb = org.Mm.eg.db,
                                  Number.of.Terms = 5,
                                  ont = "BP",
                                  pvalueCutoff = 0.05,
                                  qvalueCutoff = 0.2)
GO.exSigNet

#barplot
GO.df <- GO.exSigNet[[2]]

GO.df.1 <- data.frame(GO.term = rownames(GO.df), variable = "Gene Prop.", value = GO.df$Node.Pct, GO.name = GO.df$GO.Name)
GO.df.2 <- data.frame(GO.term = rownames(GO.df), variable = "Exp Prop.", value = GO.df$Exp.Pct, GO.name = GO.df$GO.Name)
GO.df.3 <- rbind(GO.df.1, GO.df.2)
ggplot(GO.df.3, aes(
  x = factor(GO.name, levels = unique(GO.name)),              
  y = ifelse(variable == "Gene Prop.", value, -value),  fill = variable)) +
  geom_bar(stat = 'identity', width=0.6)+  
  scale_fill_brewer(palette = 'Set2')+
  coord_flip()+                                               
  scale_y_continuous(                                         
    labels = abs,                                             
    expand = expansion(mult = c(0.1, 0.1))) +
  ylab("Proportion") +
  xlab("GO term") +
  theme_bw()

#Critical transition analysis using BioTIP-----

#load libraries
library(BioTIP)
library(cluster)
library(GenomicRanges)
library(Hmisc)
library(MASS)
require(stringr)
require(psych)
require(igraph)


meta <- data.frame(meta.data)
Data <- data.frame(Data1.AB@assays$RNA@counts)
Data <-log2(Data+1)
meta$group = meta$Type
samples <- split(meta[,"cell_barcode"],f = meta$group)
lapply(samples, length)
test <- sd_selection(Data, samples, cutoff = 0.05)

igraphL <- getNetwork(test, fdr = 0.05)

cluster <- getCluster_methods(igraphL)

names(cluster)

membersL_noweight <- getMCI(cluster, test)
plotBar_MCI(membersL_noweight,ylim = c(0,3))

maxMCIms <- getMaxMCImember(membersL_noweight[[1]],membersL_noweight[[2]], minsize = 10)
names(maxMCIms)
names(maxMCIms[[1]])
names(maxMCIms[[2]])

head(maxMCIms[['idx']])
biomodules = getMaxStats(membersL_noweight[['members']],maxMCIms[[1]])
maxMCI = getMaxStats(membersL_noweight[['MCI']],maxMCIms[[1]])
maxMCI = maxMCI[order(maxMCI,decreasing=TRUE)]
head(maxMCI)
topMCI = getTopMCI(membersL_noweight[[1]],membersL_noweight[[2]],membersL_noweight[['MCI']],min = 10)
head(topMCI)
maxSD = getMaxStats(membersL_noweight[['sd']],maxMCIms[[1]])
head(maxSD)
CTS = getCTS(topMCI, maxMCIms[[2]])
par(mar = c(10,5,0,2))
plotMaxMCI(maxMCIms,membersL_noweight[[2]],states = names(samples),las = 2)
simuMCI <- simulationMCI(length(CTS),samples,Data, B=200)
plot_MCI_Simulation(topMCI[1],simuMCI,las=2)
CTS_2 <- CTS$INCCs

IC <- getIc(Data, samples, CTS_2, PCC_sample.target = 'average')
par(mar = c(10,5,0,2))
plotIc(IC,las = 2)

# visualizing the results using a bar plot
data.IC <- data.frame(IC = IC, Group = names(IC))

data.IC %>%
  mutate(Group = fct_reorder(Group, IC, .desc = FALSE)) %>%
  ggplot(aes(x=IC, y=Group,  fill=Group)) +
    geom_bar(stat="identity", alpha=1, width=.5) +
    xlab("Index of criticality (IC)") +
    ylab("Cell group") +
    theme_bw()  

#Find critical transition signal-related signaling network-----------
#Infer the CTS-related signaling network
#setup target genes
marker.1 <- read.csv(file = "../markers/Faure2020_Mouse_marker/Faure2020_Mouse_markers_INCCs.csv")
marker.2 <- read.csv(file = "../markers/Faure2020_Mouse_marker/Faure2020_Mouse_markers_unassigned.1.csv")
marker.3 <- read.csv(file = "../markers/Faure2020_Mouse_marker/Faure2020_Mouse_markers_INCCs.BCCs.csv")

marker.1 <- marker.1$X[1:10]
marker.2 <- marker.2$X[1:10]
marker.3 <- marker.3$X[1:10]
target.123 <- c(marker.1, marker.2, marker.3)

#infer ltGRN
ltGRN <- get_ltGRN(Target = target.123, DB = exFINDER.M)

# infer potential external signals from ltGRN
Pex <- get_potentialex(Graph = ltGRN,
                Exp.Data = exp.data,
                Meta.Data = meta.data,
                cutoff = -0.05,
                AG.L = NULL)
# infer exSigNet (R, TF, T are all highly expressed, L comes from the potential external signals)
# here, we are interested in the existance of the exSigNet that associated to the CTS, so we set the cutoff values to be relatively low 
exSigNet <- get_exSigNet(Graph = ltGRN,
                       Ligands = Pex,
                       Exp.Data = exp.data,
                       Meta.Data = meta.data,
                       par.cutoff = c(0.05, 0.001), 
                       AG.R = c("INCCs", "INCCs.BCCs", "unassigned.1"),
                       AG.TF = c("INCCs", "INCCs.BCCs", "unassigned.1"))


#if overlap found, fix the overlap issue
exSigNet <- fix_OverLap(exSigNet)

# check if there is any CTS included in the inferred exSigNet
intersect(exSigNet$Graph.Node$Gene[exSigNet$Graph.Node$Role == "TF"], CTS_2)

# based on the exSigNet, we infer its sub-Network by select the CTS
Path.Graph <- get_subNet(Graph = exSigNet,
           TF = "Mitf")

# based on this "exSigNet.new", we identify the receptors that only activated by these external signals
Path.Graph <- filterLR_exSigNet(Graph = Path.Graph,
                  Exp.Data = exp.data,
                  DB = exFINDER.M)
Path.Graph <- get_NodeValue(Graph = Path.Graph,
                          Exp.Data = exp.data,
                          Meta.Data = meta.data,
                          AG.R = c("INCCs", "INCCs.BCCs", "Nociceptive"),
                          AG.TF = c("INCCs.BCCs", "Nociceptive"),
                          AG.T = c("INCCs.BCCs", "Nociceptive"),
                          Exp.ExData = NULL,
                          Meta.ExData = NULL,
                          AG.ExData = NULL)

Path.Graph <- get_EdgeWeight(Graph = Path.Graph, 
               Kh = 2)

#Visualize the network------------
library(igraph)
# customize the label for coloring
Path.Graph$Graph.Node$Role[Path.Graph$Graph.Node$Gene == "Erbb3"] = "Target (INCCs)"
Path.Graph$Graph.Node$Role[Path.Graph$Graph.Node$Gene %in% c("Rxrg", "Mpz", "Itgb3")] = "Target (INCCs.BCCs)"
Path.Graph$Graph.Node$Role[Path.Graph$Graph.Node$Gene %in% c("Myo5b", "Tgm2")] = "Target (unassigned.1)"

actors <- data.frame(
  name = Path.Graph$Graph.Node$Gene,
  role = Path.Graph$Graph.Node$Role
)

relations <- data.frame(
  from = Path.Graph$Graph.Edge$from,
  to = Path.Graph$Graph.Edge$to
)

g <- graph_from_data_frame(
  relations, directed = TRUE,
  vertices = actors
)

l <- layout_with_lgl(g) 

library(RColorBrewer)
coul  <- brewer.pal(6, "Set2") 
# Create a vector of color
my_color <- coul[as.numeric(as.factor(V(g)$role))]

plot(g, layout=l, vertex.color=my_color)
