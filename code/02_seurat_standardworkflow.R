#' ---
#' title: "scRNA-seq investigation on seaAD data"
#' author: Chiara Cabrelle 
#' ---

# Differential expression analysis results are really weird, too many DEGs so...explore the data.
# It's obvious reading the paper! They integrated data from different sources to create the atlas  
# -- Probably batch effect are the reason of our first investigation by DEA and I am using UMIs!
#

library(Seurat)
library(dplyr)
library(ggplot2)

load(file = 'results/01_seuset.rda')

# Visualize QC metrics as a violin plot
seuset[["percent.mt"]] <- PercentageFeatureSet(seuset, pattern = "^MT-")
png("plots/02_features_violin.png",w=4200,h=3200,res=300)
VlnPlot(seuset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
# ok, data were preprocessed as expected
dev.off()

# Normalize the data
seuset <- NormalizeData(seuset, normalization.method = "LogNormalize", scale.factor = 10000)
# Feature selection
seuset <- FindVariableFeatures(seuset, selection.method = "vst", nfeatures = 2000)
# Error in match.arg(arg = layer, choices = Layers(object = object, search = FALSE)) :
#   'arg' should be one of "counts", "data", "scale.data"
# no way to solve 
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seuset), 10)

# plot variable features with and without labels
png('plots/02_variablefeatures.png', w=4000,h=4000,res=300)
plot1 <- VariableFeaturePlot(seuset)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
dev.off()

# Scale the data
# all.genes <- rownames(seuset)
# seuset <- ScaleData(seuset, features = all.genes) #data are too large 
seuset <- ScaleData(seuset)
# Perform PCA
seuset <- RunPCA(seuset, features = VariableFeatures(object = seuset))
# Examine and visualize PCA results a few different ways
print(seuset[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(seuset, dims = 1:2, reduction = "pca")
png('plots/02_dimplotPCA.png', w=4000,h=4000,res=300)
DimPlot(seuset, reduction = "pca")

DimHeatmap(seuset, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(seuset, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(seuset, dims = 16:33, cells = 500, balanced = TRUE)
DimHeatmap(seuset, dims = 33:40, cells = 500, balanced = TRUE)
# Determine the dimensionality
# seuset <- JackStraw(seuset, num.replicate = 100)
# seuset <- ScoreJackStraw(seuset, dims = 1:20)
png('plots/02_elbow.png', w=3000,h=3000,res=300)
ElbowPlot(seuset, 60)
dev.off()
# cluster the cells 
seuset <- FindNeighbors(seuset, dims = 1:25)
seuset <- FindClusters(seuset, resolution = 0.8)

# Perform tsne and umap
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
seuset <- RunUMAP(seuset, dims = 1:25)
#seuset <- RunTSNE(seuset, dims = 1:30)

DimPlot(seuset, reduction = "umap")
p1 <- DimPlot(seuset, reduction = "umap", group.by = "condition")
p2 <- DimPlot(seuset, reduction = "umap", group.by = "ID_name")
png('plots/02_umap__bycond.png', w=6000,h=3000,res=300)
p1
dev.off()

png('plots/02_umap__bycelltype.png', w=3000,h=3000,res=300)
DimPlot(seuset, reduction = "umap", group.by = "celltype")
dev.off()

save(seuset, file = "results/02_seuset-DimReducted.rda")

