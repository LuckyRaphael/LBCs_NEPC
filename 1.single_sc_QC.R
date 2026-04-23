# Set library paths for Seurat version 4
.libPaths(c("C:/Users/lucky/AppData/Local/R/win-library/seurat4", .libPaths()))

# Check package versions
packageVersion("SeuratObject")
packageVersion("Seurat")

# Set system error messages to English
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())

# Set random seed
set.seed(20242025)

# Load required libraries
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(limma)
library(harmony)
library(patchwork)
library(qs)
library(data.table)

# Set working directory
setwd("G:\\2024_CRPC_NEPC\\")

# Read samples
samples = list.files('G:\\2024_CRPC_NEPC\\', pattern = 'txt')
sceList = lapply(samples, function(x) {
  print(x)
  y = file.path("G:\\2024_CRPC_NEPC\\", x)
  a = fread(y, data.table = F)
  a$Symbol <- make.unique(as.character(a$Symbol))
  a <- a[, -1]
  rownames(a) = a[, 1]
  a <- a[, -1]
  sce = CreateSeuratObject(a, project = gsub(".*_(P[0-9]+)_.*", "\\1", x), min.cells = 3, min.features = 200)
  return(sce)
})

sce.all <- merge(sceList[[1]], y = sceList[-1])
scRNA <- sce.all

scRNA@meta.data$patient <- gsub(".*_(P[0-9]+)_.*", "\\1", scRNA@meta.data$orig.ident)
scRNA@meta.data$tissue_type <- c("CRPC")

# Calculate mitochondrial and hemoglobin ratios
scRNAList <- SplitObject(scRNA, split.by = "patient")
for (i in 1:length(scRNAList)) {
  sc <- scRNAList[[i]]
  sc[["mt_percent"]] <- PercentageFeatureSet(sc, pattern = "^MT")
  HB_genes <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")
  HB_m <- match(HB_genes, rownames(sc@assays$RNA))
  HB_genes <- rownames(sc@assays$RNA)[HB_m]
  HB_genes <- HB_genes[!is.na(HB_genes)]
  sc[["HB_percent"]] <- PercentageFeatureSet(sc, features = HB_genes)
  scRNAList[[i]] <- sc
  rm(sc)
}

# Plot before filtering
violin_before <- list()
for (i in 1:length(scRNAList)) {
  violin_before[[i]] <- VlnPlot(scRNAList[[i]], group.by = "orig.ident",
                                features = c("nFeature_RNA", "nCount_RNA", "mt_percent", "HB_percent"),
                                pt.size = 0.01, ncol = 4)
}
violin_before_merge <- CombinePlots(plots = violin_before, nrow = length(scRNAList), legend = 'none', group.by = "orig.ident")
ggsave("1.violin_before_merge.pdf", plot = violin_before_merge, width = 15, height = 30)

# Merge samples and plot before filtering
scRNA <- merge(scRNAList[[1]], y = scRNAList[-1])
violin_before_scRNA <- VlnPlot(scRNA, group.by = "orig.ident",
                               features = c("nFeature_RNA", "nCount_RNA", "mt_percent", "HB_percent"),
                               pt.size = 0.01, ncol = 4)
ggsave("2.violin_before_scRNA.pdf", plot = violin_before_scRNA, width = 15, height = 7)

# Feature correlation before filtering
pdf(file = "3.featureCor.pdf", width = 22, height = 6)
plot1 <- FeatureScatter(object = scRNA, feature1 = "nCount_RNA", group.by = "orig.ident", feature2 = "mt_percent", pt.size = 0.1)
plot2 <- FeatureScatter(object = scRNA, feature1 = "nCount_RNA", group.by = "orig.ident", feature2 = "nFeature_RNA", pt.size = 0.1)
plot3 <- FeatureScatter(object = scRNA, feature1 = "mt_percent", group.by = "orig.ident", feature2 = "nFeature_RNA", pt.size = 0.1)
CombinePlots(plots = list(plot1, plot2, plot3), ncol = 3)
dev.off()

qsave(scRNA, file = "0.Merge-scRNA.qs")

# Filter cells
scRNA <- qread("G://2024_CRPC_NEPC/2.GSE137829/0.Merge-scRNA.qs")
scRNAList <- SplitObject(scRNA, split.by = "orig.ident")
scRNAList <- lapply(X = scRNAList, FUN = function(x) {
  x <- subset(x, subset = nFeature_RNA > 300 & nFeature_RNA < 7000 &
                mt_percent < 20 & HB_percent < 5 &
                nCount_RNA < quantile(nCount_RNA, 0.97) &
                nCount_RNA > 1000)
})
scRNA <- merge(scRNAList[[1]], scRNAList[2:length(scRNAList)])

# Remove doublets using scDblFinder
library(scDblFinder)
scRNA <- scRNA %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:30) %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution = 0.1)

scRNA <- as.SingleCellExperiment(scRNA)
scRNA = scDblFinder(scRNA, samples = "patient")
doublets <- DimPlot(scRNA, group.by = 'scDblFinder.class', cols = c("gold", "black", '#00bfc4', '#f8766d'))
ggsave("4.Doublets.pdf", plot = doublets, width = 10, height = 10)
scRNA <- subset(scRNA, subset = scDblFinder.class == 'singlet')

qsave(scRNA, file = "1.Doublefinder-scRNA.qs")

# DecontX for contamination correction
library(decontX)
library(scCDC)
library(Matrix)

sce.list <- SplitObject(scRNA, split.by = "orig.ident")
cor.list <- lapply(1:length(sce.list), function(i) {
  subset <- sce.list[[i]]
  if (ncol(subset@assays[["RNA"]]@data) < 1000) {
    warning("Sample has fewer than 1000 cells, skipping:", names(sce.list)[i])
    return(subset)
  }
  
  GCGs <- ContaminationDetection(subset)
  if (length(GCGs) == 0) {
    warning("No contaminated genes detected, skipping:", names(sce.list)[i])
    return(subset)
  }
  
  mislet_cont_ratio <- ContaminationQuantification(subset, rownames(GCGs))
  mislet_seuratobj_corrected <- ContaminationCorrection(subset, rownames(GCGs))
  
  mislet_seuratobj_corrected@assays[["RNA"]]@scale.data <- matrix(ncol = 0, nrow = 0)
  mislet_seuratobj_corrected@assays[["RNA"]]@data <- matrix(ncol = 0, nrow = 0)
  
  return(mislet_seuratobj_corrected)
})

scRNA <- merge(cor.list[[1]], y = cor.list[2:length(cor.list)])
scRNA@assays[["RNA_original"]] <- scRNA@assays[["RNA"]]
scRNA@assays[["RNA"]] <- scRNA@assays[["Corrected"]]
scRNA@assays[["Corrected"]] <- NULL

counts <- scRNA@assays$RNA@counts
decontX_results <- decontX(counts, z = dplyr::coalesce(scRNA$seurat_clusters, "unknown"), seed = 2024)
scRNA$contamination = decontX_results$contamination

scRNA <- scRNA %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:30) %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution = 0.1)

# Plot contamination
FeaturePlot(scRNA, features = 'contamination', raster = FALSE) +
  scale_color_viridis_c() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  xlab('UMAP_1') +
  ylab('UMAP_2')
ggsave(filename = "contamination.pdf", width = 4.2, height = 3.5)

# Filter based on contamination
scRNA <- scRNA[, scRNA$contamination < 0.2]
qsave(scRNA, file = "1.5.decontX-scRNA.qs")

# Plot violin after filtering
scRNAList <- SplitObject(scRNA, split.by = "orig.ident")
violin_scRNAter <- list()
for (i in 1:length(scRNAList)) {
  violin_scRNAter[[i]] <- VlnPlot(scRNAList[[i]], group.by = "orig.ident",
                                  features = c("nFeature_RNA", "nCount_RNA", "mt_percent", "HB_percent"),
                                  pt.size = 0.01, ncol = 4)
}
violin_scRNAter_merge <- CombinePlots(plots = violin_scRNAter, nrow = length(scRNAList), legend = 'none', group.by = "orig.ident")
ggsave("5.violin_scRNAter_merge.pdf", plot = violin_scRNAter_merge, width = 15, height = 30)

# Merge filtered samples and plot
scRNA <- merge(scRNAList[[1]], scRNAList[2:length(scRNAList)])
violin_scRNAter_scm <- VlnPlot(scRNA, group.by = "orig.ident",
                               features = c("nFeature_RNA", "nCount_RNA", "mt_percent", "HB_percent"),
                               pt.size = 0.01, ncol = 4)
ggsave("6.violin_scRNAter_scm.pdf", plot = violin_scRNAter_scm, width = 15, height = 7)

# Combine plots before and after filtering
scplot_merge <- CombinePlots(list(violin_before_scRNA, violin_scRNAter_scm), nrow = 2, legend = "none")
ggsave("7.scplot_merge.pdf", plot = scplot_merge, width = 15, height = 10)
qsave(scRNA, file = "2.QC-scRNA.qs")

# Harmony integration
scRNA <- qread("2.QC-scRNA.qs")
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 1e4)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 3000)
scRNA <- ScaleData(scRNA)
scRNA <- RunPCA(object = scRNA, pc.genes = VariableFeatures(scRNA))

# Determine dimensions for Harmony
ElbowPlot(scRNA, reduction = "pca", ndims = 50)
pct <- scRNA[["pca"]]@stdev / sum(scRNA[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
pcs <- min(which(cumu > 90 & pct < 5)[1], sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1)

# Run Harmony
library(harmony)
seuratObj <- RunHarmony(scRNA, "orig.ident")
seuratObj <- RunUMAP(seuratObj, dims = 1:pcs, reduction = "harmony")
seuratObj <- RunTSNE(seuratObj, dims = 1:pcs, reduction = "harmony")

# Visualize clusters
DimPlot(seuratObj, reduction = "umap", label = TRUE)
DimPlot(seuratObj, reduction = "tsne", label = TRUE)

scRNA <- seuratObj
scRNA <- FindNeighbors(scRNA, reduction = "harmony", dims = 1:pcs)

# Clustering with different resolutions
for (res in c(0.01, 0.05, seq(0.1, 1.5, by = 0.1))) {
  scRNA <- FindClusters(scRNA, graph.name = "RNA_snn", resolution = res)
}

# Visualize clustering results
clustree(scRNA)
ggsave(filename = "13.obj.pdf", width = 8.5, height = 10)

# Save the final object
qsave(scRNA, file = "3.harmony-scRNA.qs")

# Annotation
scRNA <- qread("3.harmony-scRNA.qs")
dir.create("./annotation")
setwd("./annotation")

# Initial annotation based on clusters
scRNA@meta.data$celltype <- NA
Idents(scRNA) <- paste0("RNA_snn_res.", resolution)

# Assign cell types based on clusters
scRNA@meta.data$celltype[scRNA@meta.data$seurat_clusters %in% c(8)] <- "Mye"
scRNA@meta.data$celltype[scRNA@meta.data$seurat_clusters %in% c(7)] <- "Mast"
scRNA@meta.data$celltype[scRNA@meta.data$seurat_clusters %in% c(18)] <- "B"
scRNA@meta.data$celltype[scRNA@meta.data$seurat_clusters %in% c(10, 11, 16)] <- "T"
scRNA@meta.data$celltype[scRNA@meta.data$seurat_clusters %in% c(19)] <- "Plasma"
scRNA@meta.data$celltype[scRNA@meta.data$seurat_clusters %in% c(3)] <- "Fib"
scRNA@meta.data$celltype[scRNA@meta.data$seurat_clusters %in% c(5)] <- "SMC"
scRNA@meta.data$celltype[scRNA@meta.data$seurat_clusters %in% c(9)] <- "Endo"
scRNA@meta.data$celltype[scRNA@meta.data$seurat_clusters %in% c(0, 1, 2, 4, 12, 15)] <- "Lum"

qsave(scRNA, file = "4.annotated_scRNA.qs")

# Visualization of markers
DefaultAssay(scRNA) <- "RNA"
Idents(scRNA) <- "celltype"

# Violin plots for selected markers
plots <- VlnPlot(scRNA, features = c("DENND4B", "PTTG1", "EPCAM", "CDH1", "KRT18", "KRT8"), 
                 group.by = "celltype", pt.size = 0, combine = FALSE)

pdf("Violin_Plots.pdf", width = 10, height = 30)
wrap_plots(plots = plots, ncol = 2)
dev.off()

# Feature plots for selected markers
featureplot <- FeaturePlot(scRNA, features = c("DENND4B", "PTTG1", "EPCAM", "CDH1", "KRT18", "KRT8"), min.cutoff = "q9")
pdf("Feature_Plots.pdf", width = 16, height = 16)
featureplot
dev.off()

# Bubble plot for selected markers
library(ggplot2)
Idents(scRNA) <- scRNA$celltype
p <- DotPlot(object = scRNA, features = features)
ggplot(p$data, aes(x = features.plot, y = id)) + 
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) + 
  facet_grid(facets = ~feature.groups, switch = "x", scales = "free_x", space = "free_x") +  
  scale_radius(breaks = c(25, 50, 75, 100), range = c(0, 4)) + 
  theme_classic() + 
  scale_color_gradient2(low = "#3C5488FF", mid = "white", high = "#DC0000FF") +
  theme(axis.text.x = element_text(angle = 90, face = 1, size = 10, hjust = 1, vjust = 0.5, color = "black"), 
        axis.text.y = element_text(size = 12, face = 1, color = "black"),
        legend.text = element_text(size = 10, face = 1), 
        legend.title = element_text(size = 15, face = 1),
        legend.position = 'top', 
        strip.placement = "outside", 
        strip.text.x = element_text(size = 15),
        axis.title = element_blank()) +
  guides(colour = guide_colourbar(title.vjust = 0.9, title.hjust = 0))

pdf("Bubble_Plot.pdf", width = 6, height = 4)
p1 + labs(tag = "Defining \ngenes of:") + 
  theme(plot.tag.position = c(0, 0.05), plot.tag = element_text(size = 12))
dev.off()