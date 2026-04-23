# =============================================================================
# PART 2: EPITHELIAL CELL SUB-ANALYSIS
# =============================================================================

setwd("F:/2025_PRAD/LYZ_NE_20250502/2.BatchRemove_BBKNN/Epi_anno/")
sce0 <- qread("merge_sce_bbknn_resolution_TME_final.qs")

# ── 1. Isolate Epithelial Cells ───────────────────────────────────────────────
sce <- sce0[, sce0$celltype == "Epi"]
DefaultAssay(sce) <- "RNA"

# ── 2. Remove Contaminating Non-Epithelial Cells (UCell Scoring) ─────────────
immune_markers <- list(
  lymphocyte = c("PTPRC","CD2","CD3D","CD3E","MS4A1","CD79A","IRF4","CD38",
                 "TNFRSF17","TPSAB1","CPA3","KIT","KLRB1","KLRC1","KLRD1"),
  myeloid    = c("LYZ","C1QA","C1QB","CSF1R","CD163","MARCO"),
  stroma     = c("VWF","SELE","PECAM1","CDH5","PLVAP","S100B")
)

counts_matrix  <- GetAssayData(sce, slot = "counts")
sparse_matrix  <- as(counts_matrix, "CsparseMatrix")
marker_score   <- ScoreSignatures_UCell(sparse_matrix, features = immune_markers)

# Retain cells with zero score across all three contaminant signatures
lymph_clean  <- rownames(marker_score)[marker_score[,"lymphocyte_UCell"] == 0]
myeloid_clean<- rownames(marker_score)[marker_score[,"myeloid_UCell"]    == 0]
stroma_clean <- rownames(marker_score)[marker_score[,"stroma_UCell"]     == 0]
pure_epi     <- Reduce(intersect, list(lymph_clean, myeloid_clean, stroma_clean))

# Additional filter: require zero expression of all contaminant markers
sce_pure <- sce[, colnames(sce) %in% pure_epi]
contam_genes  <- unique(unlist(immune_markers))
count_check   <- GetAssayData(sce_pure, assay = "RNA", layer = "counts")[
  intersect(contam_genes, rownames(sce_pure)), ]
cells_zero    <- colSums(count_check) == 0
sce <- sce_pure[, cells_zero]
qsave(sce, "Epi-sce_curated.qs")

# ── 3. Epithelial Cell Preprocessing ─────────────────────────────────────────
sce <- qread("Epi-sce_curated.qs")

# Marker gene sets for prostate epithelial subtypes
epi_markers <- list(
  AR_activity     = c("KLK3","KLK2","TMPRSS2","FKBP5","NKX3-1","PLPP1","PMEPA1","ALDH1A3","STEAP4"),
  NE_signature    = c("SYP","CHGA","CHGB","ENO2","CHRNB2","SCG3","SCN3A","PCSK1","ELAVL4","NKX2-1"),
  Squamous        = c("KRT5","KRT6A","KRT6B","DSG3","IVL","SBSN","FGFBP1","SCEL","S100A7","KRT14"),
  Basal           = c("TP63","ITGA6","KRT14","KRT15","KRT5","KRT17","CD44","CAV1","CAV2"),
  Luminal         = c("EPCAM","CDH1","TMPRSS2","NKX3-1","KLK3","KLK2","AR","KRT8","KRT18",
                      "ACPP","MSMB","AMACR","CD24","PSCA","CLDN3","ALDH1A3"),
  Club_Hillock    = c("SCGB1A1","SCGB3A1","PIGR","OLFM4","MMP7","KRT4","KRT13","LTF",
                      "SERPINB1","CLDN4","S100A16","KLF5","TACSTD2"),
  LPCs_Cycling    = c("MKI67","AURKA","DLGAP5","HMMR","TK1","CDK1","CCNB1","CCNB2",
                      "UBE2C","BIRC5","FOXM1","PLK1","NEK2"),
  NEC_markers     = c("INSM1","SOX2","ASCL1","NCAM1","CALCA","SYT1","CALCA","SCG2"),
  AR_panel        = c("AR","CHRNA2","KLK3","NKX3-1","SLC45A3","TARP"),
  Neuro_I_panel   = c("CHGA","SYP","ACTL6B","SNAP25","INSM1","ASCL1","CHRNB2","SRRM4"),
  Neuro_II_panel  = c("CELF3","PCSK1","SOX2","POU3F2","LMO3","NKX2-1"),
  squam_panel     = c("KRT5","KRT6A","KRT6B","FGFBP1"),
  EMT             = c("CDH2","FN1","VIM","TWIST1","SNAI1","ZEB1","ZEB2","DCN"),
  Stemness        = c("LTF","LY6E","TACSTD2","KRT19","PSCA","HES1","ALDH1A1","CD44","PROM1")
)

all_epi_genes <- unique(unlist(epi_markers))

sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 1e4)
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
sce <- ScaleData(sce, features = union(VariableFeatures(sce), all_epi_genes))
sce <- RunPCA(sce, features = union(VariableFeatures(sce), all_epi_genes), npcs = 50)

# PC selection
epi_var    <- sce[["pca"]]@stdev^2
epi_cumvar <- cumsum(epi_var) / sum(epi_var)
pcs        <- which(epi_cumvar >= 0.8)[1]

# ── 4. Batch Correction: Harmony + BBKNN ─────────────────────────────────────
scRNA2 <- sce %>%
  RunHarmony(group.by.vars = "patient", assay.use = "RNA", max.iter.harmony = 20) %>%
  FindNeighbors(dims = 1:pcs, reduction = "harmony") %>%
  FindClusters(resolution = 1) %>%
  RunUMAP(dims = 1:pcs, reduction = "harmony") %>%
  RunTSNE(dims = 1:pcs, reduction = "harmony")

sce_bbknn <- scRNA2 %>%
  RunBBKNN(
    reduction = "harmony", batch_key = "patient",
    n_pcs = 50, run_UMAP = TRUE, run_TSNE = TRUE,
    seed = 20252026, verbose = TRUE
  )
qsave(sce_bbknn, "Epi-sce_curated_harmony_bbknn.qs")

# ── 5. Multi-Resolution Clustering for Epithelial Cells ──────────────────────
sce <- qread("Epi-sce_curated_harmony_bbknn.qs")
sce <- FindNeighbors(sce, reduction = "harmony", dims = 1:50)

for (res in c(seq(0.01, 0.09, by = 0.01), seq(0.1, 3, by = 0.1))) {
  sce <- FindClusters(sce, graph.name = "bbknn", resolution = res, algorithm = 1)
}
qsave(sce, "Epi-sce_curated_bbknn_resolution.qs")

# ── 6. NE Scoring with NEPAL ─────────────────────────────────────────────────
NE_seurat <- NEPAL_scRNA(
  seurat.data = sce,
  method      = "AUCell",
  species     = "human",
  ncores      = 10,
  assay.names = "NE",
  DefaultAssay = FALSE
)
qsave(NE_seurat, "Epi-sce_curated_bbknn_resolution_NE.seurat.qs")

# ── 7. UCell Pathway Scoring ──────────────────────────────────────────────────
uterus <- qread("Epi-sce_curated_bbknn_resolution_NE.seurat.qs")
DefaultAssay(uterus) <- "RNA"

sparse_matrix  <- as(GetAssayData(uterus, slot = "counts"), "CsparseMatrix")
marker_score   <- ScoreSignatures_UCell(sparse_matrix, features = epi_markers)

# Add scores to metadata and as a separate assay
marker_score_df <- as.data.frame(marker_score)
rownames(marker_score_df) <- colnames(uterus)
uterus <- AddMetaData(uterus, metadata = marker_score_df)
uterus[["UCell"]] <- CreateAssayObject(data = as.matrix(t(marker_score_df)))

qsave(uterus, "Epi-sce_curated_bbknn_resolution_NE_NEPAL_Ucells.Seurat.qs")

# ── 8. Export to AnnData (.h5ad) for Python Analysis ─────────────────────────
library(reticulate); library(sceasy)
reticulate::use_python("C:/Users/lucky/anaconda3/envs/scanpy/python.exe")
use_condaenv("scanpy", required = TRUE)
loompy <- reticulate::import("loompy")

sceasy::convertFormat(
  uterus, from = "seurat", to = "anndata",
  outFile = "Epi-sce_curated_bbknn_resolution_NE_NEPAL_Ucells.h5ad"
)

# ── 9. Key Visualization ──────────────────────────────────────────────────────
FeaturePlot(uterus,
            features = c("AR","NKX3-1","TP63","SYP","CHGA","KRT5","KRT8","EPCAM"),
            order = TRUE, ncol = 4)
ggsave("Epi_key_markers.pdf", width = 16, height = 8)

ucell_cols <- grep("_UCell", colnames(uterus@meta.data), value = TRUE)
FeaturePlot(uterus, features = ucell_cols, order = TRUE, ncol = 4)
ggsave("Pathway_UCell_scores.pdf", width = 16, height = 13)