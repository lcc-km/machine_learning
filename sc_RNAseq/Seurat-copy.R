library(Seurat)
library(dplyr)
library(magrittr)
library(hdf5r)
setwd("/Users/lucc/Desktop")

pbmc <- readRDS("GSE209791_normalized_seurat_obj.rds")#该函数直接读???10x的h5文件

# store mitochondrial percentage in object meta data
pbmc<-CreateSeuratObject(count= matrix, project = "WT", min.cells = 3, min.features = 200)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")

# Show QC metrics for the first 5 cells
head(pbmc@meta.data, 5)

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

###Normalizing the data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

###Identification of highly variable features (feature selection)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10)
plot1 + plot2
plot2
###Scaling the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

###Perform linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
# PCA
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

###Determine the ‘dimensionality??? of the dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)

###Cluster the cells
pbmc <- FindNeighbors(pbmc, dims = 1:15)
pbmc <- FindClusters(pbmc, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)






###Run non-linear dimensional reduction (UMAP/tSNE)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:14)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
p1 <- DimPlot(pbmc, reduction = "umap",label = TRUE)
p1
#saveRDS(pbmc, file = "../pbmc_tutorial.rds")

pbmc <- RunTSNE(pbmc, dims = 1:15)
head(pbmc@reductions$tsne@cell.embeddings)
p2 <- DimPlot(pbmc, reduction = "tsne",label = TRUE)
p2
#saveRDS(pbmc, file = "pbmc_tutorial.rds") 

### Finding differentially expressed features (cluster biomarkers)
# find all markers of cluster 1
cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 10)

cluster3.markers <- FindMarkers(pbmc, ident.1 = 3, min.pct = 0.25)
head(cluster3.markers, n = 10)

cluster4.markers <- FindMarkers(pbmc, ident.1 = 4, min.pct = 0.25)
head(cluster4.markers, n = 10)

cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, min.pct = 0.25)
head(cluster5.markers, n = 10)

cluster6.markers <- FindMarkers(pbmc, ident.1 = 6, min.pct = 0.25)
head(cluster6.markers, n = 10)

cluster7.markers <- FindMarkers(pbmc, ident.1 = 7, min.pct = 0.25)
head(cluster7.markers, n = 10)

cluster8.markers <- FindMarkers(pbmc, ident.1 = 8, min.pct = 0.25)
head(cluster8.markers, n = 10)

cluster9.markers <- FindMarkers(pbmc, ident.1 = 9, min.pct = 0.25)
head(cluster9.markers, n = 10)

cluster10.markers <- FindMarkers(pbmc, ident.1 = 10, min.pct = 0.25)
head(cluster10.markers, n = 10)

cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, min.pct = 0.25)
head(cluster0.markers, n = 10)




# find all markers
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25)
top3 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
DoHeatmap(pbmc,features = top3$gene) + NoLegend()

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# run sctransform
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)

# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

write.table(pbmc.markers,"pbmc.markers_1",sep = "\t")
write.table((pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)),"pbmc.markers_2",sep = "\t")

cluster1.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(pbmc, features = c("tekt3", "lfng", "sost", "fndc7rs4", "atoh1b",
                           "dld", "ebf3a", "isl1","ovgp1","sfrp1a",
                           "si:ch73-261i21.5","wnt11r"), 
        pt.size = 0.2, ncol = 4, slot = "counts", log = TRUE)

FeaturePlot(pbmc, features = c("tekt3", "lfng", "sost", "fndc7rs4", "atoh1b",
                               "dld", "ebf3a", "isl1","ovgp1","sfrp1a",
                               "si:ch73-261i21.5","wnt11r"),
            pt.size = 0.2, ncol = 3)

###Hair Cell lineage
FeaturePlot(pbmc, features = c("tekt3","atoh1b","dld"),
            pt.size = 0.2, ncol = 3)
###Central Cells
FeaturePlot(pbmc, features = c("lfng","ebf3a","isl1"),
            pt.size = 0.2, ncol = 3)

### Polar Cells
FeaturePlot(pbmc, features = c("sost"),
            pt.size = 0.2, ncol = 3)

###Mantle Cells
FeaturePlot(pbmc, features = c("ovgp1","sfrp1a"),
            pt.size = 0.2, ncol = 3)

###Ring (not in HC lineage)
FeaturePlot(pbmc, features = c("fndc7rs4"),
            pt.size = 0.2, ncol = 3)

###A-P Cells
FeaturePlot(pbmc, features = c("si:ch73-261i21.5","wnt11r"),
            pt.size = 0.2, ncol = 3)


top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()


### Assigning cell type identity to clusters #####不知道没有做
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()



















# These are now standard steps in the Seurat workflow for visualization and clustering
pbmc <- RunPCA(pbmc, verbose = FALSE)
pbmc <- RunUMAP(pbmc, dims = 1:30, verbose = FALSE)

pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
pbmc <- FindClusters(pbmc, verbose = FALSE)
DimPlot(pbmc, label = TRUE) + NoLegend()

VlnPlot(pbmc, features = c("tekt3", "lfng", "sost", "fndc7rs4", "atoh1b",
                           "did", "ebf3a", "isl1","ovgp1","sfrp1a",
                           "si:ch73-261i21.5","wnt11r"), 
        pt.size = 0.2, ncol = 4)
FeaturePlot(pbmc, features = c("tekt3", "lfng", "sost", "fndc7rs4", "atoh1b",
                               "did", "ebf3a", "isl1","ovgp1","sfrp1a",
                               "si:ch73-261i21.5","wnt11r"),
            pt.size = 0.2, ncol = 3)
