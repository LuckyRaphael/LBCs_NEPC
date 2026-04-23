# =============================================================================
# inferCNV Analysis Pipeline: CNV Detection in Prostate Cancer Epithelial Cells
# =============================================================================

# ── 0. Environment Setup ─────────────────────────────────────────────────────
.libPaths(c("C:/Users/lucky/AppData/Local/R/win-library/seurat4", .libPaths()))
rm(list = ls())
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 80 * 1024 * 1024^2)
set.seed(20242025)

library(Seurat); library(qs); library(infercnv); library(rjags)
library(ComplexHeatmap); library(circlize); library(RColorBrewer)
library(ggplot2); library(ggpubr); library(ggalluvial); library(ggbeeswarm)
library(dplyr); library(tidyverse); library(patchwork)

# =============================================================================
# STEP 1: Prepare Input Files for inferCNV
# =============================================================================

# ── 1.1 Load Epithelial Seurat Object & Reference Cells ──────────────────────
sce_CNV <- qread("Epi-sce_curated_bbknn_resolution_NE_NEPAL_Ucells_Filtered_final.Seurat.qs")
load("reference_mat.Rdata")   # loads spike_mat1 (Endo) and spike_mat2 (Mast)

# ── 1.2 Build Combined Count Matrix (Epithelial + Reference cells) ────────────
epiMat    <- as.data.frame(GetAssayData(sce_CNV, slot = "counts", assay = "RNA"))
shared    <- intersect(rownames(epiMat), rownames(spike_mat1))
this_dat  <- cbind(epiMat[shared, ], spike_mat1[shared, ], spike_mat2[shared, ])
save(this_dat, file = "sce_CNV_counts.matrix.rda")

# ── 1.3 Build Cell Annotation File ───────────────────────────────────────────
# Use resolution 1.5 clusters as epithelial group labels
Idents(sce_CNV)      <- sce_CNV$bbknn_res.1.5
sce_CNV$epi_cluster  <- paste0("Epi_", as.character(Idents(sce_CNV)))

groupinfo <- data.frame(
  cell_id = colnames(this_dat),
  group   = c(
    sce_CNV$epi_cluster,          # epithelial tumor cells
    rep("spike-1", ncol(spike_mat1[, 1:2000])),   # Endo spike
    rep("ref-1",   ncol(spike_mat1[, 2001:3000])),# Endo reference
    rep("spike-2", ncol(spike_mat2[, 1:2000])),   # Mast spike
    rep("ref-2",   ncol(spike_mat2[, 2001:3000])) # Mast reference
  )
)
write.table(groupinfo, file = "groupFiles.txt",
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# =============================================================================
# STEP 2: Run inferCNV
# =============================================================================

load("sce_CNV_counts.matrix.rda")

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = this_dat,
  annotations_file  = "groupFiles.txt",
  delim             = "\t",
  gene_order_file   = "hg38_gencode_v27.txt",   # hg38 gene position file
  ref_group_names   = c("ref-1", "ref-2"),       # normal reference cells
  chr_exclude       = c("chrX", "chrY", "chrM")  # exclude sex/mito chromosomes
)

# Notes on ref_group selection:
#   1. Best: matched normal tissue samples
#   2. Alternative: immune cells (if confirmed non-malignant)
#   3. Last resort: omit ref_group (uses all-cell mean as baseline)

infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff            = 0.1,          # 0.1 for 10x Genomics; 1 for Smart-seq2
  out_dir           = "Epi-Infercnv",
  no_prelim_plot    = TRUE,
  cluster_by_groups = TRUE,
  denoise           = TRUE,
  HMM               = FALSE,
  min_cells_per_gene = 10,
  num_threads       = 18,
  write_expr_matrix = TRUE          # exports infercnv.observations.txt
)
save(infercnv_obj, file = "infercnv_obj.rda")

# ── Optional: Custom heatmap with modified colors ─────────────────────────────
load("infercnv_obj.rda")
infercnv::plot_cnv(
  infercnv_obj,
  output_filename  = "inferCNV_heatmap",
  output_format    = "pdf",
  custom_color_pal = color.palette(c("#2067AE", "white", "#B11F2B"))
)

# =============================================================================
# STEP 3: CNV Score Quantification
# =============================================================================

# ── 3.1 Method 1 — KMeans Clustering on CNV Expression Matrix ────────────────
# Reference: inferCNV heatmap visual inspection — cells with high CNV = malignant

load("sce_CNV_counts.matrix.rda")
load("infercnv_obj.rda")

expr   <- infercnv_obj@expr.data
expr1  <- t(expr)
save(expr1, file = "expr1.rda")

# Relabel reference groups with descriptive names
groupinfo         <- read.table("groupFiles.txt")
cell_anno         <- groupinfo
colnames(cell_anno) <- c("cell_id", "group")
cell_anno$group[cell_anno$group == "spike-1"] <- "spike_Endo"
cell_anno$group[cell_anno$group == "ref-1"]   <- "ref_Endo"
cell_anno$group[cell_anno$group == "spike-2"] <- "spike_Mast"
cell_anno$group[cell_anno$group == "ref-2"]   <- "ref_Mast"

# KMeans clustering (k = 5; adjust based on data)
set.seed(2025)
load("expr1.rda")
kmeans_result <- kmeans(expr1, centers = 5)
kmeans_df     <- data.frame(
  cell_id   = rownames(expr1),
  k_cluster = as.factor(kmeans_result$cluster)
) %>%
  inner_join(cell_anno, by = "cell_id") %>%
  arrange(k_cluster)

# Row annotation for heatmap
annotation_row <- data.frame(
  k_cluster = kmeans_df$k_cluster,
  group     = kmeans_df$group,
  row.names = kmeans_df$cell_id
)

# Cluster composition summary
result <- annotation_row %>%
  group_by(k_cluster, group) %>%
  summarise(count = n(), .groups = "drop")
print(result)

# ── 3.2 Custom Color Palette ──────────────────────────────────────────────────
color_palette <- c(
  "#96cb8f", "#8bc96d", "#4dae47", "#5c9e43", "#a4cde1",
  "#67a4cc", "#277fb8", "#549da3", "#b79973", "#f38989",
  "#ec5051", "#e32427", "#ef6a45", "#f9b769", "#f9a341",
  "#f48521", "#ee8e46", "#d4a6a8", "#af93c4", "#8660a8",
  "#815e99", "#c6b598", "#f6f28f", "#d4a55b", "#b05a28"
)

epi_colors <- c(
  "Epi_0"  = "#A4CDE1", "Epi_1"  = "#277FB8", "Epi_2"  = "#4DAE47",
  "Epi_3"  = "#EC5051", "Epi_4"  = "#B79973", "Epi_5"  = "#F6F28F",
  "Epi_6"  = "#EE8E46", "Epi_7"  = "#EF6A45", "Epi_8"  = "#96CB8F",
  "Epi_9"  = "#815E99", "Epi_10" = "#F9A341", "Epi_11" = "#AF93C4",
  "Epi_12" = "#B05A28"
)

# ── 3.3 ComplexHeatmap with Chromosome Annotation ────────────────────────────
gene_pos   <- read.delim("hg38_gencode_v27.txt", header = FALSE)
gene_pos   <- gene_pos[gene_pos$V1 %in% rownames(expr), ]
new_cluster <- unique(gene_pos$V2)

top_color <- HeatmapAnnotation(
  cluster = anno_block(
    labels    = gsub("chr", "", new_cluster),
    gp        = gpar(col = "white"),
    labels_gp = gpar(cex = 1, col = "black"),
    height    = unit(5, "mm")
  )
)

n_clusters    <- length(unique(kmeans_df$k_cluster))
color_cluster <- colorRampPalette(ggsci::pal_npg()(10))(n_clusters)
names(color_cluster) <- as.character(seq_len(n_clusters))

left_anno <- rowAnnotation(
  df  = annotation_row,
  col = list(
    group     = c(epi_colors,
                  "ref_Endo"   = "grey20", "spike_Endo" = "grey40",
                  "spike_Mast" = "grey60", "ref_Mast"   = "grey80"),
    k_cluster = color_cluster
  ),
  show_annotation_name = FALSE
)

# ── 3.4 Method 2 — CNV Score (sum of squared deviations from 1) ──────────────
# Reference: Peng et al., Cell Research (pancreatic PDAC heterogeneity)
CNV_score <- as.data.frame(colMeans((expr - 1)^2))
colnames(CNV_score) <- "CNV_score"
CNV_score$cell_id   <- rownames(CNV_score)
CNV_score           <- CNV_score %>% inner_join(kmeans_df, by = "cell_id")
save(CNV_score, file = "CNV_score.rda")

# Add CNV score to Seurat metadata
sce_CNV <- AddMetaData(sce_CNV, CNV_score)

# Violin plot: CNV score by epithelial cluster
VlnPlot(sce_CNV, features = "CNV_score", group.by = "epi_cluster",
        pt.size = 0) +
  scale_fill_manual(values = epi_colors) +
  labs(x = NULL, y = "CNV Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# CNV score violin by k-means cluster
p_kmeans <- ggplot(CNV_score, aes(k_cluster, log(CNV_score))) +
  geom_violin(aes(fill = k_cluster), color = "black") +
  scale_fill_manual(values = color_cluster) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1),
               geom = "pointrange", color = "black", size = 0.5) +
  theme_bw() +
  theme(axis.text   = element_text(colour = "black", size = 12),
        axis.title  = element_text(color = "black", size = 12),
        panel.grid  = element_blank())
ggsave("CNV_score_kmeans.pdf", p_kmeans, width = 5, height = 4)

# CNV score violin by epithelial group
p_group <- ggplot(CNV_score, aes(group, log(CNV_score))) +
  geom_violin(aes(fill = group), color = "black") +
  scale_fill_manual(values = epi_colors) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1),
               geom = "pointrange", color = "black", size = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title  = element_text(color = "black", size = 12),
        panel.grid  = element_blank())

p_kmeans + p_group
ggsave("CNV_score_cluster_group.pdf", width = 9.2, height = 3.6)

# Sankey diagram: k-means cluster vs epithelial group
summary_data <- CNV_score %>%
  group_by(k_cluster, group) %>%
  summarise(count = n(), .groups = "drop")

ggplot(summary_data, aes(axis1 = k_cluster, axis2 = group, y = count)) +
  geom_alluvium(aes(fill = group)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  labs(x = "K-means Cluster", y = "Cell Count")
ggsave("CNV_Sankey.pdf", width = 8, height = 8)

# =============================================================================
# STEP 4: CNV Score + Correlation Method (Malignant Cell Identification)
# =============================================================================
# Reference: Moncada et al., Nature Biotechnology (pancreatic cancer TME)

system.time({
  obs <- read.table("./Epi-Infercnv/infercnv.observations.txt", header = TRUE, check.names = FALSE)
  ref <- read.table("./Epi-Infercnv/infercnv.references.txt",   header = TRUE, check.names = FALSE)
})

estimateCNV <- function(obs, ref, score_threshold, cor_threshold) {
  # CNV score: mean squared deviation from copy-neutral baseline (1)
  cnv_obs <- colMeans((obs - 1)^2)
  cnv_ref <- colMeans((ref - 1)^2)

  # Select top 5% highest-CNV cells as reference signature
  cell_top <- names(sort(cnv_obs, decreasing = TRUE))[1:round(length(cnv_obs) * 0.05)]
  cnv_top  <- rowMeans(obs[, cell_top])

  # Pearson correlation of each cell vs. top-CNV signature
  cor_obs <- apply(obs, 2, function(x) cor(x, cnv_top))
  cor_ref <- apply(ref, 2, function(x) cor(x, cnv_top))

  cnv <- data.frame(
    score   = c(cnv_obs, cnv_ref),
    cor     = c(cor_obs, cor_ref),
    cell_id = c(colnames(obs), colnames(ref))
  )

  # Classify cells
  cnv$type <- "Other"
  cnv$type[cnv$score >  score_threshold & cnv$cor >  cor_threshold] <- "Malignant"
  cnv$type[cnv$score <= score_threshold & cnv$cor <= cor_threshold] <- "Not Malignant"
  return(cnv)
}

CNV_score_cor <- estimateCNV(obs, ref, score_threshold = 0.001, cor_threshold = 0.25)
colnames(CNV_score_cor)[3] <- "cell_id"
CNV_score_cor <- CNV_score_cor %>% inner_join(kmeans_df, by = "cell_id")
save(CNV_score_cor, file = "CNV_score_cor_EEC.rda")

# =============================================================================
# STEP 5: Data-Driven Malignancy Threshold (Based on Benign Cell Distribution)
# =============================================================================

load("CNV_score_cor_EEC.rda")

# Use benign cells to set one-tailed 95% CI thresholds
benign_df <- CNV_score_cor %>% filter(group %in% c("ref_Endo", "ref_Mast"))

confidence_level <- 0.95
fil_score <- quantile(benign_df$score, confidence_level)
fil_cor   <- quantile(benign_df$cor,   confidence_level)
cat("CNV Score threshold:", fil_score, "\nCNV Cor threshold:", fil_cor, "\n")

# Reclassify all cells using benign-derived thresholds
CNV_score_cor <- CNV_score_cor %>%
  mutate(type = case_when(
    score > fil_score & cor > fil_cor ~ "Malignant",
    score < fil_score & cor < fil_cor ~ "Not Malignant",
    TRUE ~ "Other"
  ))

# Add to Seurat metadata
sce_CNV <- AddMetaData(sce_CNV, CNV_score_cor %>%
                         select(cell_id, score, cor, type) %>%
                         column_to_rownames("cell_id"))

# =============================================================================
# STEP 6: Visualization
# =============================================================================

# ── 6.1 CNV Score vs Correlation Scatter Plots ───────────────────────────────
load("CNV_score_cor_EEC.rda")

p_by_group <- ggplot(CNV_score_cor, aes(score, cor, color = group)) +
  geom_point(size = 0.05) +
  scale_color_manual(values = epi_colors) +
  labs(x = "CNV Score", y = "CNV Correlation") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(1, 0), legend.justification = c("right", "bottom"),
        legend.background = element_rect(color = "black")) +
  guides(color = guide_legend(override.aes = list(size = 3)))

p_by_type <- ggplot(CNV_score_cor, aes(score, cor, color = type)) +
  geom_point(size = 0.05) +
  scale_color_manual(values = c("Malignant" = "#009E73", "Not Malignant" = "black", "Other" = "grey")) +
  labs(x = "CNV Score", y = "CNV Correlation") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(1, 0), legend.justification = c("right", "bottom"),
        legend.background = element_rect(color = "black")) +
  guides(color = guide_legend(override.aes = list(size = 3)))

p_by_group + p_by_type
ggsave("CNV_score_cor_scatter.pdf", width = 10, height = 5)

# ── 6.2 UMAP Colored by Malignancy Prediction ────────────────────────────────
malignant_colors <- c("Malignant" = "#FFD300", "Not Malignant" = "black", "Other" = "grey50")

df_umap <- sce_CNV@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  mutate(
    cluster    = sce_CNV$type,
    epi_cluster = sce_CNV$epi_cluster,
    tissue     = sce_CNV$Tissue,
    cnv_score  = sce_CNV$score,
    cnv_cor    = sce_CNV$cor
  )

p_umap <- ggplot(df_umap, aes(UMAP_1, UMAP_2, color = cluster)) +
  geom_point(size = 0.01, alpha = 0.8, shape = 16) +
  scale_color_manual(values = malignant_colors) +
  facet_wrap(~ tissue, scales = "free", ncol = 4) +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        legend.title = element_blank(),
        legend.text  = element_text(size = 14)) +
  guides(fill = guide_legend(override.aes = list(size = 5)))
ggsave("Epi_CNV_UMAP_malignant.pdf", p_umap, width = 6.2, height = 3)

# ── 6.3 UMAP Colored by CNV Score (Continuous Gradient) ──────────────────────
p_score <- ggplot(df_umap %>% arrange(cnv_score),
                  aes(UMAP_1, UMAP_2, colour = cnv_score)) +
  geom_point(size = 0.01, alpha = 0.8, shape = 16) +
  scale_color_gradientn(
    colors = c("black", "#333300", "#DAA520", "#FFD300", "gold"),
    values = scales::rescale(c(2.806e-05, 150e-05, 350e-05, 600e-05, 1026e-05))
  ) +
  theme_classic() +
  theme(axis.text  = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 15))
ggsave("Epi_CNV_UMAP_score.pdf", p_score, width = 3.8, height = 3)

# ── 6.4 Beeswarm Plot: CNV Score Ordered by Median per Cluster ───────────────
# Order clusters by median CNV score
cluster_order <- df_umap %>%
  group_by(epi_cluster) %>%
  summarise(med = median(cnv_score, na.rm = TRUE)) %>%
  arrange(med) %>%
  pull(epi_cluster)

df_tumor <- df_umap %>%
  filter(tissue == "Tumor") %>%
  mutate(
    epi_cluster = factor(epi_cluster, levels = cluster_order),
    highlight   = ifelse(cluster == "Malignant", "Malignant", as.character(epi_cluster))
  )

p_beeswarm <- ggplot(df_tumor, aes(epi_cluster, cnv_score)) +
  geom_quasirandom(aes(color = highlight), size = 0.0001, alpha = 0.7) +
  scale_color_manual(values = c("Malignant" = "black", epi_colors)) +
  labs(x = "Epithelial Cluster", y = "CNV Score") +
  theme_classic() +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y  = element_text(size = 12),
        axis.title   = element_text(size = 15),
        legend.position = "none")
ggsave("CNV_beeswarm_malignant.pdf", p_beeswarm, width = 8, height = 4)