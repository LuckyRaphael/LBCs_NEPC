# =============================================================================
# scRNA-seq Analysis Pipeline: Batch Correction, Clustering & Epithelial Annotation
# Project: PRAD (Prostate Adenocarcinoma) NE Subtype Analysis
# =============================================================================

# ── 0. Environment Setup ─────────────────────────────────────────────────────
Sys.setenv(LANGUAGE = "en")          # English error messages
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 80 * 1024 * 1024^2)  # 80 GB memory limit
set.seed(2024)
rm(list = ls())

# Load Seurat v4 from dedicated library path (coexist with v5)
.libPaths(c("C:/Users/lucky/AppData/Local/R/win-library/seurat4", .libPaths()))

library(Seurat); library(qs); library(dplyr); library(patchwork)
library(harmony); library(bbknnR); library(UCell); library(NEPAL)
library(ggplot2); library(ggsci); library(reshape2)

# =============================================================================
# PART 1: WHOLE-TUMOR MICROENVIRONMENT (TME) ANALYSIS
# =============================================================================

# ── 1. Load Merged & Filtered Data ───────────────────────────────────────────
setwd("F:/2025_PRAD/LYZ_NE_20250502/2.BatchRemove_BBKNN/")
system.time({
  scRNA <- qread("F:/2025_PRAD/LYZ_NE_20250502/1.Merge_Filtered_Metadata/scRNA_merge(r-raw)_filtered_gene_metadata.qs")
})

# ── 2. Normalization & Variable Feature Selection ─────────────────────────────
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 1e4)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
scRNA <- ScaleData(scRNA)
scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA))

# ── 3. PCA Diagnostics ───────────────────────────────────────────────────────
pdf("8.PCAheatmap.pdf", width = 20, height = 30)
DimHeatmap(scRNA, dims = 1:30, cells = 500, balanced = TRUE, nfeatures = 30, ncol = 3)
dev.off()

pdf("9.PCAGene.pdf", width = 10, height = 30)
VizDimLoadings(scRNA, dims = 1:30, reduction = "pca", nfeatures = 30)
dev.off()

# ── 4. Determine Optimal PC Number ───────────────────────────────────────────
# Method 1: Elbow + cumulative variance threshold
ElbowPlot(scRNA, reduction = "pca", ndims = 50)
ggsave("10.1.Elbow.pdf", width = 10, height = 10)

pct  <- scRNA[["pca"]]@stdev / sum(scRNA[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1  <- which(cumu > 90 & pct < 5)[1]           # >90% cumulative variance, <5% per PC
co2  <- sort(which((pct[-length(pct)] - pct[-1]) > 0.1), decreasing = TRUE)[1] + 1

# Method 2: Variance ratio
variances  <- scRNA@reductions$pca@stdev^2
cum_var    <- cumsum(variances) / sum(variances)
n_pcs      <- which(cum_var >= 0.8)[1]
print(paste("Suggested PCs:", n_pcs))

pcs <- n_pcs   # use 80% variance threshold

# Elbow plot with threshold annotation
plot_df <- data.frame(pct = pct, cumu = cumu, rank = seq_along(pct))
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) +
  geom_text() + geom_vline(xintercept = 80, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") + theme_bw()
ggsave("10.2.Elbow.pdf", width = 6, height = 5)

# ── 5. Batch Correction: BBKNN ───────────────────────────────────────────────
sce_bbknn <- scRNA %>%
  RunBBKNN(
    reduction = "pca",
    batch_key  = "patient",
    n_pcs      = 50,
    run_UMAP   = TRUE,
    seed       = 2025,
    verbose    = TRUE
  )
qsave(sce_bbknn, "merge_sce_bbknn_resolution_CRPC_NEPC.qs")

# ── 6. Multi-Resolution Clustering ───────────────────────────────────────────
scRNA <- qread("merge_sce_bbknn_resolution_CRPC_NEPC.qs")

# Recalculate variance-based PC cutoff
variances <- scRNA@reductions$pca@stdev^2
n_pcs     <- which(cumsum(variances) / sum(variances) >= 0.8)[1]
pcs       <- n_pcs

sce <- scRNA
sce <- FindNeighbors(sce, reduction = "pca", dims = 1:pcs)

for (res in c(seq(0.01, 0.09, by = 0.01), seq(0.1, 3, by = 0.1))) {
  sce <- FindClusters(sce, graph.name = "bbknn", resolution = res, algorithm = 2)
  # Algorithm: 1=Louvain, 2=Louvain+multilevel, 3=SLM, 4=Leiden
}

# Visualize resolution sweep (UMAP)
res_list <- c(seq(0.01, 0.09, by = 0.01), seq(0.1, 3, by = 0.1))
plot_panels <- lapply(res_list, function(r) {
  DimPlot(sce, reduction = "umap",
          group.by = paste0("bbknn_res.", r), label = TRUE) & NoAxes()
})
cluster_umap <- wrap_plots(plotlist = plot_panels, ncol = 4)
ggsave("7.umap_resolution.pdf", cluster_umap, width = 20, height = 35)
qsave(sce, "merge_sce_bbknn_resolution.qs")

# ── 7. Cell Type Marker Dictionary ───────────────────────────────────────────
CellMarkers <- list(
  "T_NK"    = c("CD3D","CD4","IL7R","CD8A","CD8B","GZMB","PRF1","GZMA","IFNG",
                "FOXP3","IL2RA","CTLA4","NKG7","KLRK1","CD56","CD16"),
  "B_Plasma"= c("CD19","MS4A1","IGHD","IGHM","TCL1A","CD27","IGHG1","IGHA1",
                "SDC1","CD38","PRDM1","XBP1","MZB1","JCHAIN"),
  "CAF"     = c("DCN","PDGFRA","PDPN","ACTA2","POSTN","FAP","MMP11","CTHRC1",
                "IL6","CFD","CXCL12","CD74","HLA-DR"),
  "Pericyte"= c("PDGFRB","RGS5","KCNJ8","CD248"),
  "Macro"   = c("CD68","CD163","MRC1","CSF1R","CD86","IL10","ARG1","CD14"),
  "Mono"    = c("CD14","FCGR3A","LYZ","S100A12","CCR2","VCAN"),
  "Neut"    = c("CSF3R","MPO","CXCR2","S100A8","S100A9","CEACAM8"),
  "Endo"    = c("PECAM1","CD34","VWF","CDH5","KDR","CLDN5","LYVE1","FLT4"),
  "Mast"    = c("KIT","FCER1A","TPSAB1","CPA3","HDC","MS4A2"),
  "DCs"     = c("CLEC9A","XCR1","CD1C","CLEC10A","LILRA4","IRF8","HLA-DRA"),
  "Epi_Lum" = c("KRT8","KRT18","EPCAM","NKX3-1","KLK2","KLK3","AR","ACPP","MSMB"),
  "Epi_Bas" = c("TP63","KRT5","KRT14","KRT15"),
  "Epi_NE"  = c("SYP","CHGA","CHGB","ENO2","NCAM1","SCG2","CALCA")
)

# Generate DotPlots across all resolutions for each cell type
setwd("F:/2025_PRAD/LYZ_NE_20250502/2.BatchRemove_BBKNN/cell-anno-resolution/")
for (cell in names(CellMarkers)) {
  markers <- unique(unlist(CellMarkers[[cell]]))
  pdf(paste0(cell, "_Resolution_all_dotplot.pdf"), width = 18, height = 8)
  for (res in res_list) {
    Idents(sce) <- sce@meta.data[[paste0("bbknn_res.", res)]]
    p <- DotPlot(sce, features = markers) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      labs(subtitle = paste0("resolution = ", res))
    print(p)
  }
  dev.off()
}

# ── 8. Manual Cell Type Annotation (resolution = 3) ──────────────────────────
resolution <- 3
Idents(sce) <- sce@meta.data[[paste0("bbknn_res.", resolution)]]
sce@meta.data$celltype <- NA

# Assign cell types by cluster identity
sce@meta.data$celltype[Idents(sce) %in% c(34)]         <- "NK cell"
sce@meta.data$celltype[Idents(sce) %in% c(3, 13, 18)]  <- "CTLs"
sce@meta.data$celltype[Idents(sce) %in% c(10)]         <- "CD8+ T cell"
sce@meta.data$celltype[Idents(sce) %in% c(1)]          <- "CD4+ T cell"
sce@meta.data$celltype[Idents(sce) %in% c(5)]          <- "Tregs"
sce@meta.data$celltype[Idents(sce) %in% c(21)]         <- "B cell"
sce@meta.data$celltype[Idents(sce) %in% c(30)]         <- "Plasma"
sce@meta.data$celltype[Idents(sce) %in% c(9)]          <- "Pericyte"
sce@meta.data$celltype[Idents(sce) %in% c(4)]          <- "myCAF"
sce@meta.data$celltype[Idents(sce) %in% c(23)]         <- "iCAF"
sce@meta.data$celltype[Idents(sce) %in% c(8, 20)]      <- "Mono/Macro"
sce@meta.data$celltype[Idents(sce) %in% c(33)]         <- "DCs"
sce@meta.data$celltype[Idents(sce) %in% c(26)]         <- "Neut"
sce@meta.data$celltype[Idents(sce) %in% c(7, 15)]      <- "Endo"
sce@meta.data$celltype[Idents(sce) %in% c(24)]         <- "Mast"

# Epithelial cells annotated at lower resolution for broader capture
Idents(sce) <- sce@meta.data[["bbknn_res.0.08"]]
sce@meta.data$celltype[Idents(sce) %in% c(0)] <- "Epi"

# Label remaining unassigned cells as "Other" and remove
sce@meta.data$celltype[is.na(sce@meta.data$celltype)] <- "Other"
sce <- subset(sce, celltype != "Other")

DimPlot(sce, group.by = "celltype", label = TRUE)
setwd("F:/2025_PRAD/LYZ_NE_20250502/2.BatchRemove_BBKNN/Epi_anno/")
ggsave("TME_annotation.pdf", width = 5.5, height = 4)

Idents(sce) <- sce$celltype
qsave(sce, "merge_sce_bbknn_resolution_TME_final.qs")

# ── 9. TME DotPlot (Publication Figure) ──────────────────────────────────────
curated_features <- list(
  "Lymphoid" = c("CD3D","CD4","IL7R","CD8A","CD8B","GZMA","IFNG",
                 "NKG7","KLRK1","CTLA4","IL2RA","FOXP3",
                 "MS4A1","CD19","TCL1A","JCHAIN","TNFRSF17","CD38"),
  "Stromal"  = c("COL1A2","DCN","CTHRC1","POSTN","FAP","PDGFRA","ACTA2",
                 "MYH11","IL6","RGS5","KCNJ8"),
  "Myeloid"  = c("CD163","CD86","CD68","CD14","MRC1","CSF1R","LYZ","CCR2",
                 "IL10","S100A12","CSF3R","IL8","CD1C","CLEC10A","FCER1A"),
  "Endothel" = c("PECAM1","CD34","VWF","CDH5","KDR","CLDN5"),
  "Mast"     = c("KIT","TPSAB1","CPA3","HDC","MS4A2","TPSB2"),
  "Epithelial"= c("KRT8","KRT18","EPCAM","NKX3-1","KLK2","KLK3",
                  "AR","ACPP","MSMB","KRT5","KRT14","KRT15")
)

uterus <- sce
Idents(uterus) <- uterus$celltype
levels(uterus) <- c("CD4+ T cell","CD8+ T cell","CTLs","NK cell","Tregs",
                    "B cell","Plasma",
                    "myCAF","iCAF","Pericyte",
                    "Mono/Macro","Neut","DCs",
                    "Endo","Mast","Epi")

# Styled DotPlot
p <- DotPlot(uterus, features = curated_features) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

p_styled <- ggplot(p$data, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  facet_grid(~ feature.groups, switch = "x", scales = "free_x", space = "free_x") +
  scale_radius(breaks = c(25, 50, 75, 100), range = c(0, 4)) +
  scale_color_gradient2(low = "#3C5488FF", mid = "white", high = "#DC0000FF") +
  theme_classic() +
  theme(
    axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10, color = "black"),
    axis.text.y  = element_text(size = 12, color = "black"),
    strip.placement = "outside",
    strip.text.x = element_text(size = 15),
    legend.position = "top",
    axis.title   = element_blank()
  ) +
  labs(size = "Percent Expressed", color = "Average Expression")

pdf("TME_DotPlot_styled.pdf", width = 12, height = 5)
print(p_styled)
dev.off()