library("ArchR")
library(BSgenome.Mmusculus.UCSC.mm10)
library("parallel")

rm(list = ls())

setwd("./result_E13_neg2")
genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.Mmusculus.UCSC.mm10)
genomeAnnotation

inputFiles <- c(
                "./Input_data/GSM5024885_E13_MGE_dlxneg_2_ATAC_fragments.tsv.gz"
)

names(inputFiles) <- c("GSM5024885_E13_MGE_dlxneg_2")


addArchRThreads(threads = 10)
addArchRGenome("mm10")

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4,
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
ArrowFiles
class(ArrowFiles)

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, 
  knnMethod = "UMAP", 
  LSIMethod = 1
)

projHeme1 <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "HemeTutorial",
  copyArrows = TRUE 
)

projHeme1
paste0("Memory Size = ", round(object.size(projHeme1) / 10^6, 3), " MB")
getAvailableMatrices(projHeme1)

head(projHeme1$cellNames)
head(projHeme1$Sample)
quantile(projHeme1$TSSEnrichment)


projHeme1[1:100, ]

projHeme1[projHeme1$cellNames[1:100], ]

# sample name
idxSample <- BiocGenerics::which(projHeme1$Sample %in% "GSM5024885_E13_MGE_dlxneg_2")
cellsSample <- projHeme1$cellNames[idxSample]
projHeme1[cellsSample, ]
# TSS enrichment score
idxPass <- which(projHeme1$TSSEnrichment >= 8)
cellsPass <- projHeme1$cellNames[idxPass]
projHeme1[cellsPass, ]

df <- getCellColData(projHeme1, select = "nFrags")
df <- getCellColData(projHeme1, select = c("log10(nFrags)", "nFrags - 1"))

df <- getCellColData(projHeme1, select = c("log10(nFrags)", "TSSEnrichment"))
p <- ggPoint(
  x = df[,1],
  y = df[,2],
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
  ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")

p

p1 <- plotGroups(
  ArchRProj = projHeme1,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "TSSEnrichment",
  plotAs = "ridges"
)
p1
plot(1:5,1:5) 

p2 <- plotGroups(
  ArchRProj = projHeme1,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
p2

p3 <- plotGroups(
  ArchRProj = projHeme1,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "log10(nFrags)",
  plotAs = "ridges"
)
p3

p4 <- plotGroups(
  ArchRProj = projHeme1,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
p4

plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics_E13_neg2.pdf", ArchRProj = projHeme1, addDOC = FALSE, width = 4, height = 4)

p1 <- plotFragmentSizes(ArchRProj = projHeme1)
p1

p2 <- plotTSSEnrichment(ArchRProj = projHeme1)
p2

plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile_E13_neg2.pdf", ArchRProj = projHeme1, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = projHeme1, outputDirectory = "Save-ProjHeme1", load = FALSE)

projHeme2 <- filterDoublets(projHeme1)
projHeme2

projHemeTmp <- filterDoublets(projHeme1, filterRatio = 1.5)

rm(projHemeTmp)

projHeme2 <- addIterativeLSI(
  ArchRProj = projHeme2,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30
)

projHeme2 <- addClusters(
  input = projHeme2,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8
)
head(projHeme2$Clusters)
table(projHeme2$Clusters)

cM <- confusionMatrix(paste0(projHeme2$Clusters), paste0(projHeme2$Sample))
cM

library(pheatmap)
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
  mat = as.matrix(cM), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black"
)
p

projHeme2 <- addClusters(
  input = projHeme2,
  reducedDims = "IterativeLSI",
  method = "scran",
  name = "ScranClusters",
  k = 15
)

projHeme2 <- addUMAP(
  ArchRProj = projHeme2, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)

p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters_E13_neg2.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)
ggsave("./result_E13_neg2/HemeTutorial/Plots/Plot-UMAP1-Sample-Clusters_E13_neg2.png", p1, width = 5, height = 5, dpi = 300)
ggsave("./result_E13_neg2/HemeTutorial/Plots/Plot-UMAP2-Sample-Clusters_E13_neg2.png", p2, width = 5, height = 5, dpi = 300)

projHeme2 <- addTSNE(
  ArchRProj = projHeme2, 
  reducedDims = "IterativeLSI", 
  name = "TSNE", 
  perplexity = 30
)
p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "TSNE")
p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "TSNE")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-TSNE-Sample-Clusters_E13_neg2.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)
ggsave("./result_E13_neg2/HemeTutorial/Plots/Plot-TSNE1-Sample-Clusters_E13_neg2.png", p1, width = 5, height = 5, dpi = 300)
ggsave("./result_E13_neg2/HemeTutorial/Plots/Plot-TSNE2-Sample-Clusters_E13_neg2.png", p2, width = 5, height = 5, dpi = 300)

markersGS <- getMarkerFeatures(
  ArchRProj = projHeme2, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
markerList$C6

###可视化
markerGenes  <- c(
  "Egr4"
)

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = markerGenes,
  transpose = TRUE
)

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = markerGenes,
  transpose = TRUE
)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap_E13_neg2", width = 8, height = 6, ArchRProj = projHeme2, addDOC = FALSE)


markerGenes  <- c(
  "Egr4"
)

p <- plotEmbedding(
  ArchRProj = projHeme2, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAP",
  quantCut = c(0.01, 0.95),
  imputeWeights = NULL
)
p$Tac1
p2 <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))

plotPDF(plotList = p, 
        name = "UMAP_Marke_E13_neg2_Plot-UMAP-Marker-Genes-WO-Imputation.pdf", 
        ArchRProj = projHeme2, 
        addDOC = FALSE, width = 5, height = 5)
ggsave("./result_E13_neg2/HemeTutorial/Plots/UMAP_Marke_E13_neg2_Plot-UMAP-Marker-Genes-WO-Imputation.png", p, width = 5, height = 5, dpi = 300)

projHeme3 = projHeme2
projHeme4 <- addGroupCoverages(ArchRProj = projHeme3, groupBy = "Clusters")

pathToMacs2 = "/home/miniconda3/bin/macs2"

projHeme4 <- addReproduciblePeakSet(
  ArchRProj = projHeme4, 
  groupBy = "Clusters",
  pathToMacs2 = "/home/miniconda3/bin/macs2"
)
getPeakSet(projHeme4)

saveArchRProject(ArchRProj = projHeme4, outputDirectory = "Save-ProjHeme4", load = FALSE)

projHeme5 <- addPeakMatrix(projHeme4)
getAvailableMatrices(projHeme5)


t <- getMatrixFromProject(projHeme5,useMatrix ="GeneScoreMatrix" )
colData(t)
gene <- rowData(t)
metadata(t)
tt <- as.data.frame(assay(t))
tt <- as.matrix(assay(t))
rownames(tt) <- gene$name
write.csv(tt,"GSM5024885_E13_MGE_dlxneg_2.GeneScoreMatrix.csv",row.names = T)

pt <- getMatrixFromProject(projHeme5,useMatrix ="PeakMatrix" )
peak <- rowRanges(pt)
metadata(pt)
ptt <- as.data.frame(assay(pt))
ptt <- as.matrix(assay(pt))
rownames(ptt) <- peak$idx
#peak <- ranges(rowranges)
peak <- as.data.frame(peak)
ptt <- cbind(peak,ptt)
write.csv(ptt,"GSM5024885_E13_MGE_dlxneg_2.PeakMatrix.csv",row.names = T)
