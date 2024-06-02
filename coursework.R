library(Seurat)
library(SeuratData)
library(patchwork)
library(Matrix)
library(patchwork)
library(ggplot2)
library(magrittr)
library(dplyr)
library(celldex)
library(SingleR)
library(pheatmap)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SingleR")
BiocManager::install("celldex")


data <- Read10X(data.dir = "/Users/annaorlova/Desktop/scRNAcoursework/GSE185231_RAW")
data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

sdata <- CreateSeuratObject(counts = data)

# Calculating the percentage of mitochondrial genes
sdata[["percent.mt"]] <- PercentageFeatureSet(sdata, pattern = "^MT-")

####QC####
VlnPlot(sdata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(sdata, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sdata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

plot1 <- plot1 + scale_y_continuous(limits = c(0, 25))

plot2 <- plot2 + scale_y_continuous(limits = c(0, 10000))

plot1 + plot2
sdata <- subset(sdata, subset = nFeature_RNA < 7500 & percent.mt < 15)




sdata <- NormalizeData(sdata)
sdata <- FindVariableFeatures(sdata, selection.method = "vst", nfeatures = 7500)


top10 <- head(VariableFeatures(sdata), 10)
top30 <- head(VariableFeatures(sdata), 30)
top30

# plot variable features with and without labels
plot3 <- VariableFeaturePlot(sdata)
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE)
plot3 + plot4


all.genes <- rownames(sdata)
sdata <- ScaleData(sdata, features = all.genes)

sdata <- RunPCA(sdata, features = VariableFeatures(object = sdata))


####Clusterization#####
#print(sdata[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(sdata, dims = 1:2, reduction = "pca")
DimPlot(sdata, reduction = "pca")
DimHeatmap(sdata, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(sdata)

sdata <- FindNeighbors(sdata, dims = 1:13)
sdata <- FindClusters(sdata, resolution = c(0.03, 0.05, 0.07, 0.1,0.2, 0.3, 0.5, 0.7, 1))
View(sdata@meta.data)
DimPlot(sdata, group.by = "RNA_snn_res.0.03", label = TRUE)
sdata <- FindClusters(sdata, resolution = 0.03)

sdata <- RunUMAP(sdata, dims = 1:13)
DimPlot(sdata, reduction = "umap")

#####defining cell types####
ref <- celldex::HumanPrimaryCellAtlasData()
View(as.data.frame(colData(ref)))


sdata_counts <- GetAssayData(sdata, slot = 'counts')

pred <- SingleR(test = sdata_counts,
                ref = ref,
                labels = ref$label.main)

pred

sdata$singleR.labels <- pred$labels[match(rownames(sdata@meta.data), rownames(pred))]
DimPlot(sdata, reduction = 'umap', group.by = 'singleR.labels')

pred
pred$scores

View(as.data.frame(pred$scores))

plotScoreHeatmap(pred)

plotDeltaDistribution(pred)

tab <- table(Assigned=pred$labels, Clusters=sdata$seurat_clusters)
pheatmap(log10(tab+10), color = colorRampPalette(c('white','blue'))(10))



###from paper
gene_exp_plt8 <- FeaturePlot(sdata, features = c("CD3D", "CD68", "PTPRC", "MZB1", "MOG", "IDH1"), reduction = 'umap', 
                             max.cutoff = 3, ncol = 3)
gene_exp_plt8

gene_exp_plt4 <- FeaturePlot(sdata, features = c("PDGFRA", "GFAP", "CDK4", 'NES', 'CCND1', 'MDM2'), reduction = 'umap', 
                             max.cutoff = 3, ncol = 3)
gene_exp_plt4

###from NCG High-grade glioma

gene_exp_plt <- FeaturePlot(sdata, features = c("ACVR1", "TP53", "KRAS",  "PTEN", "ATRX", "PIK3CA"), reduction = 'umap', 
                            max.cutoff = 3, ncol = 3)
gene_exp_plt

gene_exp_plt3 <- FeaturePlot(sdata, features = c("NF1", "PIK3R1", "PPM1D", "BCOR", "BCORL1", "BRAF"), reduction = 'umap', 
                             max.cutoff = 3, ncol = 3)
gene_exp_plt3

gene_exp_plt2 <- FeaturePlot(sdata, features = c("CCND2", "CCND3", "CDK6",  "CDKN2A", "EGFR", "MET"), reduction = 'umap', 
                            max.cutoff = 3, ncol = 3)
gene_exp_plt2

gene_exp_plt5 <- FeaturePlot(sdata, features = c("MYC", "MYCN", "NTRK1",  "NTRK3", "KMT2D", "CDKN2C"), reduction = 'umap', 
                             max.cutoff = 3, ncol = 3)
gene_exp_plt5

DefaultAssay(sdata) <- 'prediction.score.celltype'
#########Markers DE####
markers_all = FindAllMarkers(sdata,genes.use = VariableFeatures(sdata),
                             only.pos = TRUE, 
                             min.pct = 0.25, 
                             thresh.use = 0.25)

markers_all_single <- markers_all[markers_all$gene %in% names(table(markers_all$gene))[table(markers_all$gene) == 1],]

table(table(markers_all_single$gene))
table(markers_all_single$cluster)
head(markers_all_single)

cluster0.markers <- FindMarkers(sdata, ident.1 = 0)
head(cluster0.markers, n = 30)

cluster1.markers <- FindMarkers(sdata, ident.1 = 1)
head(cluster1.markers, n = 30)


markers_all %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 15) %>%
  ungroup() -> top15
DoHeatmap(sdata, features = top15$gene) + NoLegend()

###by clusters###
###Macrophages###

Macrophages <- FeaturePlot(sdata, features = c("TREM2", "TNF", "CD14",  "CD163", "CD40", "CD68"), reduction = 'umap', 
                             max.cutoff = 3, ncol = 3)
Macrophages

###T-cells###

tcells <- FeaturePlot(sdata,  features = c("CTLA4", "CD5", "IL5", "CD27", "CD7", "STAT4"), reduction = 'umap',
                      max.cutoff = 3, ncol = 3)

tcells


###B-cells - Same as t-cells###

###Oligodendricitoes###


new.cluster.ids <- c("Tumor cells", "Macrophages", "T-cells", "Oligodendrocytes")
# Переименовываем кластеры
names(new.cluster.ids) <- levels(sdata)
sdata_2 <- RenameIdents(sdata, new.cluster.ids)

# Проверяем изменения
Idents(sdata_2)
umap_plot <- DimPlot(sdata_2, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("UMAP Plot with Renamed Clusters") +
  theme_minimal()

# Отображаем график
print(umap_plot)



VlnPlot(sdata, features = c("TP53", "IDH1", 'PDGFRA', 'NF1', 'KRAS', 'PIK3CA'))

markers_all %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
DoHeatmap(sdata, features = top20$gene) + NoLegend()



#####POLARIZATION OF MACROPHAGES####

macro_cells <- WhichCells(sdata, idents = 1)
macro_subset <- subset(sdata, cells = macro_cells)




VizDimLoadings(macro_subset, dims = 1:2, reduction = "pca")
DimPlot(macro_subset, reduction = "pca")
DimHeatmap(macro_subset, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(macro_subset)

macrophages <- FindNeighbors(macro_subset, dims = 1:9)
macrophages <- FindClusters(macro_subset, resolution = c(0.03, 0.05, 0.07, 0.1,0.2, 0.3, 0.5, 0.7, 1))
#View(macrophages@meta.data)
DimPlot(macrophages, group.by = "RNA_snn_res.0.1", label = TRUE)
macrophages <- FindClusters(macrophages, resolution = 0.1)

macrophages <- RunUMAP(macrophages, dims = 1:13)
DimPlot(macrophages, reduction = "umap")

macro_markers_all <- FindAllMarkers(macrophages, only.pos = TRUE)
macro_markers_all %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 15) %>%
  ungroup() -> top15
DoHeatmap(macrophages, features = top15$gene) + NoLegend()

M1 <- FeaturePlot(macrophages, features = c("CD86", "CD68", "CD80", "TLR2", "TLR4", "SOCS3"), reduction = 'umap', 
                             max.cutoff = 3, ncol = 3)
M1

M1_ <- FeaturePlot(macrophages, features = c("TNF", "CCL3", "CCL4", "CCL5", "CCL2", "MHS"), reduction = 'umap', 
                    max.cutoff = 3, ncol = 3)
M1_

M2 <- FeaturePlot(macrophages, features = c("CD4", "MRC1", "STAT6", "HGF", "TLR8", "TLR1"), reduction = 'umap', 
                    max.cutoff = 3, ncol = 3)
M2

new.macro.ids <- c("M1-Macrophages", "M2-Macrophages", "M3-Macrophages")

names(new.macro.ids) <- levels(macrophages)
macrophages <- RenameIdents(macrophages, new.macro.ids)

Idents(macrophages)
umap_plot_macro <- DimPlot(macrophages, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("UMAP Plot with Renamed Clusters of Macrophages") +
  theme_minimal()


print(umap_plot_macro)

VlnPlot(macrophages, features = c("TP53", "CD44", 'PDGFRA', 'CCND1', 'MIF', 'SPP1'))


macro_clusters <- Idents(macrophages)
macro_clusters <- data.frame(Cells = names(macro_clusters), Cluster = as.character(macro_clusters))


sdata$new_macro_cluster <- NA


for (cell in macro_clusters$Cells) {
  sdata$new_macro_cluster[cell] <- macro_clusters$Cluster[macro_clusters$Cells == cell]
}


sdata$new_macro_cluster <- ifelse(is.na(sdata$new_macro_cluster), as.character(Idents(sdata)), sdata$new_macro_cluster)


Idents(sdata) <- sdata$new_macro_cluster


DimPlot(sdata, reduction = "umap", group.by = "new_macro_cluster", label = TRUE)


DimPlot(sdata, reduction = "umap", label = FALSE)



####GLIOMA CELLS SUBPOPULATION#####

cancer_cells <- WhichCells(sdata, idents = 0)
cancer_subset <- subset(sdata, cells = cancer_cells)

ElbowPlot(cancer_subset)
cancer <- FindNeighbors(cancer_subset, dims = 1:15)
cancer <- FindClusters(cancer_subset, resolution = c(0.03, 0.05, 0.07, 0.1,0.2, 0.25, 0.15))
#View(cancer@meta.data)
DimPlot(cancer , group.by = "RNA_snn_res.0.2", label = TRUE)
cancer <- FindClusters(cancer, resolution = 0.2)

cancer <- RunUMAP(cancer, dims = 1:15)
DimPlot(cancer, reduction = "umap")

cancer1 <- FeaturePlot(cancer, features = c("TP53", "OLIG1", "BCAN", "PLP1", "SOX4", "DLL3" ), reduction = 'umap',
                              max.cutoff = 3, ncol = 3)
cancer1

cancer2 <- FeaturePlot(cancer, features = c("PDGFRA", "OLIG1", "BCAN", "PLP1", "SOX4", "DLL3" ), reduction = 'umap',
                       max.cutoff = 3, ncol = 3)
cancer2

cancer_markers_all <- FindAllMarkers(cancer, only.pos = TRUE)
cancer_markers_all %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 15) %>%
  ungroup() -> top15
DoHeatmap(cancer, features = top15$gene) + NoLegend()

new.cancer.ids <- c("OPC", "Glial Support", "Active Proliferation", 
                   "High Survival and Inflammation", "Migration and Adhesion", "Immune Inflammatory Response")

names(new.cancer.ids) <- levels(cancer)
cancer <- RenameIdents(cancer, new.cancer.ids)


Idents(cancer)
umap_plot_cancer <- DimPlot(cancer, reduction = "umap", label = FALSE, pt.size = 0.5) +
  ggtitle("UMAP Plot with Renamed Clusters of cancer") +
  theme_minimal()


print(umap_plot_cancer)

VlnPlot(cancer, features = c("TP53", "CD44", 'PDGFRA', 'CCND1', 'MIF', 'SPP1'))

VlnPlot(cancer, features = c("GIMAP7", "BTG1", 'TNFAIP3', 'PTEN', 'MIF', 'SPP1'))
VlnPlot(macrophages, features = c("GIMAP7", "BTG1", 'TNFAIP3', 'PTEN', 'MIF', 'SPP1'))
VlnPlot(immune, features = c("GIMAP7", "BTG1", 'TNFAIP3', 'PTEN', 'MIF', 'SPP1'))

FeaturePlot(sdata, features = c("BTG1", 'TNFAIP3', 'PTEN', 'MIF'), reduction = 'umap',
            max.cutoff = 3, ncol = 3)

FeaturePlot(sdata, features = c("PDGFRA", 'NF1', 'MET', 'CDK4', 'KRAS', 'GFAP'), reduction = 'umap',
            max.cutoff = 3, ncol = 3)




cancer_clusters <- Idents(cancer)
cancer_clusters <- data.frame(Cells = names(cancer_clusters), Cluster = as.character(cancer_clusters))


sdata$new_cancer_cluster <- NA


for (cell in cancer_clusters$Cells) {
  sdata$new_cancer_cluster[cell] <- cancer_clusters$Cluster[cancer_clusters$Cells == cell]
}


sdata$new_cancer_cluster <- ifelse(is.na(sdata$new_cancer_cluster), as.character(Idents(sdata)), sdata$new_cancer_cluster)


Idents(sdata) <- sdata$new_cancer_cluster


DimPlot(sdata, reduction = "umap", group.by = "new_cancer_cluster", label = TRUE)


DimPlot(sdata, reduction = "umap", label = FALSE)



#####T-CELLS VS B-CELLS#####
immune_cells <- WhichCells(sdata, idents = 2)
immune_subset <- subset(sdata, cells = immune_cells)

ElbowPlot(immune_subset)
immune <- FindNeighbors(immune_subset, dims = 1:15)
immune <- FindClusters(immune_subset, resolution = c(0.03, 0.05, 0.07, 0.1,0.2, 0.25, 0.3))
#View(cancer@meta.data)
DimPlot(immune , group.by = "RNA_snn_res.0.2", label = TRUE)
immune <- FindClusters(immune, resolution = 0.2)

immune <- RunUMAP(immune, dims = 1:15)
DimPlot(immune, reduction = "umap")

tcells <- FeaturePlot(immune, features = c("CD2", "CD3", "CD4", "CD5", "CD7", "CD8"), reduction = 'umap', 
                  max.cutoff = 3, ncol = 3)
tcells

bcells <- FeaturePlot(immune, features = c("CD19", "CD20", "CD22", "IgM", "CD27", "CD38"), reduction = 'umap', 
                   max.cutoff = 3, ncol = 3)
bcells

immune_map<- FeaturePlot(immune, features = c("BTG1", "GIMAP4", "GIMAP7", "TNFAIP3"), reduction = 'umap', 
                         max.cutoff = 3, ncol = 3)

immune_map

immune_markers_all <- FindAllMarkers(immune, only.pos = TRUE)
immune_markers_all %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 15) %>%
  ungroup() -> top15
DoHeatmap(immune, features = top15$gene) + NoLegend()

new.immune.ids <- c("B-cells", "T-cells", "Plasma cells")

names(new.immune.ids) <- levels(immune)
immune <- RenameIdents(immune, new.immune.ids)


Idents(immune)
umap_plot_immune <- DimPlot(immune, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("UMAP Plot with Renamed Clusters of cancer") +
  theme_minimal()


print(umap_plot_immune)



immune_clusters <- Idents(immune)
immune_clusters <- data.frame(Cells = names(immune_clusters), Cluster = as.character(immune_clusters))


sdata$new_cluster <- NA


for (cell in immune_clusters$Cells) {
  sdata$new_cluster[cell] <- immune_clusters$Cluster[immune_clusters$Cells == cell]
}


sdata$new_cluster <- ifelse(is.na(sdata$new_cluster), as.character(Idents(sdata)), sdata$new_cluster)


Idents(sdata) <- sdata$new_cluster


DimPlot(sdata, reduction = "umap", group.by = "new_cluster", label = TRUE)


sdata <- RenameIdents(sdata, `3` ="Oligodendrocytes")
DimPlot(sdata, reduction = "umap", label = TRUE)
DimPlot(sdata, reduction = "umap", label = FALSE)

#####LIGAND-RECEPTOR CELL INTERACTION#####
install.packages("circlize")
install.packages("Rcpp")
install.packages("RcppArmadillo")

install.packages("/Users/annaorlova/Desktop/scRNAcoursework/CellChat-main", repos = NULL, type = "source")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")
BiocManager::install("NMF")

devtools::install_local("/Users/annaorlova/Desktop/scRNAcoursework/CellChat-main")
library(CellChat)


cellchat <- createCellChat(object = sdata, group.by = "ident")


CellChatDB <- CellChatDB.human

cellchat@DB <- CellChatDB



cellchat <- subsetData(cellchat)


devtools::install_github('immunogenomics/presto')


cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)

cellchat <- identifyOverExpressedInteractions(cellchat)


cellchat <- projectData(cellchat, PPI.human)


cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)


df.net <- subsetCommunication(cellchat)



cellchat <- computeCommunProbPathway(cellchat)

# visualizing
cellchat <- aggregateNet(cellchat)

netVisual_circle(cellchat@net$count, vertex.weight = as.numeric(table(cellchat@idents)), weight.scale = T, label.edge= F, title.name = "Number of interactions")

netVisual_heatmap(cellchat)

netAnalysis_contribution(cellchat, signaling = "TNF")

# avalible signaling checking 
available_signaling <- unique(cellchat@netP$pathway)
print(available_signaling)

netAnalysis_contribution(cellchat, signaling = "CCL")
netAnalysis_contribution(cellchat, signaling = "SPP1")
netVisual_heatmap(cellchat, signaling = "SPP1", color.heatmap = "Reds", title.name = "SPP1 signaling network")
netVisual_heatmap(cellchat, signaling = "CCL", color.heatmap = "Reds", title.name = "CCL signaling network")

netVisual_heatmap(cellchat, signaling = "APP", color.heatmap = "Reds", title.name = "APP signaling network")




netVisual_heatmap(cellchat, signaling = "MIF", color.heatmap = "Reds", title.name = "MIF signaling network")




cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
#bubble plots for signal networks 

tumor_cells <-  c("OPC", "Glial Support", "Active Proliferation", 
                                                    "High Survival and Inflammation", "Migration and Adhesion", 
                                                    "Immune Inflammatory Response")
stromal_cells <- c("M1-Macrophages", "M2-Macrophages", "M3-Macrophages", 
                                    "Oligodendrocytes", "T-cells", "B-cells", "Plasma cells")
macros <-  c("M1-Macrophages", "M2-Macrophages", "M3-Macrophages")

immunes <-  c( "Oligodendrocytes", "T-cells", "B-cells", "Plasma cells")
netVisual_bubble(cellchat, 
                 sources.use = macros,  
                 targets.use = tumor_cells,   
                 signaling = NULL,   
                 remove.isolate = TRUE,
                 color.heatmap = "Spectral") 

netVisual_bubble(cellchat, 
                 sources.use = immunes,  
                 targets.use = tumor_cells,   
                 signaling = NULL,   
                 remove.isolate = TRUE,
                 color.heatmap = "Spectral") 

netVisual_bubble(cellchat, 
                 sources.use = "Oligodendrocytes",  
                 targets.use = tumor_cells,   
                 signaling = NULL,   
                 remove.isolate = TRUE,
                 color.heatmap = "Spectral") 

netVisual_bubble(cellchat, 
                 sources.use = stromal_cells,   
                 targets.use = tumor_cells,  
                 signaling = "ApoE",     
                 remove.isolate = FALSE, 
                 color.heatmap = "Spectral") 

netVisual_bubble(cellchat, 
                 sources.use = tumor_cells,   
                 targets.use = stromal_cells,  
                 signaling = "JAM",     
                 remove.isolate = FALSE, 
                 color.heatmap = "Spectral") 

netVisual_bubble(cellchat, 
                 sources.use = tumor_cells,   
                 targets.use = stromal_cells,  
                 signaling = "NRXN",     
                 remove.isolate = FALSE, 
                 color.heatmap = "Spectral") 


netVisual_bubble(cellchat, 
                 sources.use = stromal_cells,   
                 targets.use = tumor_cells,   
                 signaling = "SPP1",     
                 remove.isolate = FALSE, 
                 color.heatmap = "Spectral") 

netVisual_bubble(cellchat, 
                 sources.use = stromal_cells,   
                 targets.use = tumor_cells,   
                 signaling = "MIF",     
                 remove.isolate = FALSE, 
                 color.heatmap = "Spectral") 
netVisual_bubble(cellchat, 
                 sources.use = stromal_cells,   
                 targets.use = tumor_cells,   
                 signaling = "CD99",     
                 remove.isolate = FALSE, 
                 color.heatmap = "Spectral") 

netVisual_heatmap(cellchat, signaling = "ApoE", color.heatmap = "Reds", title.name = "APOE signaling network")
