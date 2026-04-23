# ################# 1. Seurat 5 Installation #######################
.libPaths()
packageVersion("Seurat")

# ################# 2. Seurat 4 Installation #######################
.libPaths(c("C:/Users/lucky/AppData/Local/R/win-library/seurat4", .libPaths()))
packageVersion("Seurat")

# ################# 3. Error Handling ########################
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)

# ################# 4. Set Random Seed ########################
set.seed(2025)

# ################# 5. Multi-threading Setup ######################
options(future.globals.maxSize = 70 * 1024^2)  # Set global size limit to 70GB
memory.limit(size=70000)  # Set memory limit to 70GB

# ################# 6. Load Required Libraries ###################
library(Seurat)
library(qs)

# ################# Epithelial Cell Annotation ###################
setwd("F://2025_PRAD/LYZ_NE_20250502/2.BatchRemove_BBKNN/Epi_anno/Epi_reanno/hdWGCNA/")
sce <- qread("Epi-sce_curated_bbknn_resolution_NE_NEPAL_Ucells_Filtered_final.Seurat.qs")
sce$celltype <- sce$celltype2
sce$sample <- sce$patient

DimPlot(sce, group.by = "celltype")

# Extract target cell names
Target_cells <- rownames(sce@meta.data)[
  !sce@meta.data$celltype %in% c("Basal cell", "DNPC", "ARLPC_1", "ARPC_1") & 
  !sce@meta.data$tissue_type %in% c("Treatment-naive")
]

# Create a Seurat object with target cells
Target_seurat <- subset(sce, cells = Target_cells)
DefaultAssay(Target_seurat) <- "RNA"

# Count cells per patient
patient_cell_counts <- table(Target_seurat$patient)
patients_to_keep <- names(patient_cell_counts[patient_cell_counts >= 100])
Target_seurat <- Target_seurat[, Target_seurat$patient %in% patients_to_keep]

DimPlot(Target_seurat, group.by = "celltype")

umap_coords <- Target_seurat@reductions$umap@cell.embeddings
filtered_cells <- Target_seurat@meta.data[umap_coords[, 1] >= -1 & umap_coords[, 2] >= -1, ]
Target_seurat <- Target_seurat[, !colnames(Target_seurat) %in% row.names(filtered_cells)]
DimPlot(Target_seurat)

# #################### hdWGCNA Setup ############################
setwd("F://2025_PRAD/LYZ_NE_20250502/2.BatchRemove_BBKNN/Epi_anno/Epi_reanno/hdWGCNA/")
library(hdWGCNA)
library(WGCNA)
library(tidyverse)

set.seed(2025520)
enableWGCNAThreads(nThreads = 8)

sce <- Target_seurat
DefaultAssay(sce) <- 'RNA'

# ################# 1. Pre-processing ############################
sce <- SetupForWGCNA(sce,
                     wgcna_name = "Bra",
                     gene_select = "fraction",
                     fraction = 0.01)

# ################# 2. Construct Metacells ######################
sce <- MetacellsByGroups(
  sce,
  group.by = c("celltype", "orig.ident"),
  k = 20,
  reduction = 'umap',
  ident.group = 'celltype'
)

sce <- NormalizeMetacells(sce)

# ################# 3. Co-expression Network Analysis ###########
sce <- SetDatExpr(sce, assay = 'RNA', slot = 'data')
sce <- TestSoftPowers(sce, networkType = 'signed')

# Plot soft powers
plot <- PlotSoftPowers(sce)
wrap_plots(plot, ncol = 2)
ggsave("1.soft-power_threshold.pdf", width = 7, height = 5)

# Construct co-expression network
sce <- ConstructNetwork(
  overwrite_tom = TRUE,
  sce,
  soft_power = 9,
  tom_name = "NEPC",
  setDatExpr = FALSE,
  corType = "pearson",
  networkType = "signed",
  TOMType = "signed",
  tom_outdir = "TOM"
)

# ################# Module Identification #######################
pdf(file = "2.scRNA_hdWGCNA_Dendrogram.pdf", width = 4, height = 2.5)
PlotDendrogram(sce, main = 'scRNA hdWGCNA Dendrogram')
dev.off()

# Compute module eigengenes
sce <- ScaleData(sce)
sce <- ModuleEigengenes(sce, scale.model.use = "linear", group.by.vars = 'patient')

# Compute module connectivity
sce <- ModuleConnectivity(sce, corFnc = "bicor", corOptions = "use='p'", harmonized = TRUE)

# Get the module table
modules <- GetModules(sce)
write.csv(modules, file = 'modules.csv')

# ################# Visualization ###############################
# Plot genes ranked by kME for each module
PlotKMEs(sce, ncol = 3)

# Plot module eigengenes
plot_hMEs <- ModuleFeaturePlot(sce, reduction = "umap", features = 'MEs', order = TRUE, raster = TRUE)
pdf("5_ModuleFeaturePlot_MEs.pdf", width = 10, height = 8)
wrap_plots(plot_hMEs, ncol = 5)
dev.off()

# Save the Seurat object
save(sce, file = 'sce_wgcna.RData')