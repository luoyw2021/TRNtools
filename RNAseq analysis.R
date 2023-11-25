library("ArchR")
library(BSgenome.Mmusculus.UCSC.mm10)
library("parallel")
library('Seurat')
library(dplyr)
library(patchwork)

rm(list = ls())

#dir="./input_data/GSM5014302_MGE_Dlx6pos/"
#dir="./input_data/GSM5014303_MGE_Dlx6neg"
#dir="./input_data/GSM5014304_P2_Cortex_Dlxpos"
#dir="./input_data/GSM5014305_P10_ALM_Dlxpos"
dir="./input_data/GSM5014307_P28_ALM_Dlxpos"

pbmc.counts <- Read10X(data.dir = dir)

pbmc <- CreateSeuratObject(counts = pbmc.counts)
pbmc
str(pbmc)

pbmc <- CreateSeuratObject(counts = pbmc.counts, project = "pbmc3k", min.cells = 3, min.features = 200)


pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(pbmc), 10)

plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1+plot2

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)


pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)


pbmc <- FindNeighbors(pbmc, dims = 1:10)

pbmc <- FindClusters(pbmc, resolution = 0.5)

head(Idents(pbmc), 5)


pbmc <- RunUMAP(pbmc, dims = 1:10, label = T)
head(pbmc@reductions$umap@cell.embeddings) 

p1 <- DimPlot(pbmc, reduction = "umap")

pbmc <- RunTSNE(pbmc, dims = 1:10)
head(pbmc@reductions$tsne@cell.embeddings)
p2 <- DimPlot(pbmc, reduction = "tsne")
p1 + p2
saveRDS(pbmc, file = "pbmc_tutorial.rds")  


cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25)
?FindMarkers
top3 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
DoHeatmap(pbmc, features = top3$gene) + NoLegend()


pFeaturePlot <- FeaturePlot(pbmc, features = c("Egr4", "Gad2", "Lhx6", "Tac1", "Pvalb", "Sst"), pt.size = 0.25)
ggsave("/output/FeaturePlot_Marker_P28_ALM.png", pFeaturePlot, dpi = 300)

pVlnPlot <- VlnPlot(pbmc, features = c("Egr4", "Gad2", "Lhx6", "Tac1", "Pvalb", "Sst"), slot = "counts", log = TRUE)
ggsave("/output/VlnPlot_Marker_P28_ALM.png", pVlnPlot, dpi = 300)

