# ==============================================================================
# Single-cell RNA-seq Data Preprocessing Pipeline
# ==============================================================================
# Description: Comprehensive preprocessing pipeline for 10X Genomics scRNA-seq data
# Author: Jason
# Date: 7/23/2025
# ==============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

# Set global parameters
set.seed(12)
RANDOM_SEED <- 12
PROJECT_NAME <- "PBMC_scRNAseq"
DATA_DIR <- "data/filtered_feature_bc_matrix"
FIGURES_DIR <- "figures"
RESULTS_DIR <- "results"
DPI <- 300

# Create output directories if they don't exist
if (!dir.exists(FIGURES_DIR)) {
  dir.create(FIGURES_DIR, recursive = TRUE)
}
if (!dir.exists(RESULTS_DIR)) {
  dir.create(RESULTS_DIR, recursive = TRUE)
}

# ==============================================================================
# 1. DATA LOADING AND INITIAL SETUP
# ==============================================================================

message("Loading 10X Genomics data...")
raw_counts <- Read10X(data.dir = DATA_DIR)
message(paste("Loaded", nrow(raw_counts), "genes and", ncol(raw_counts), "cells"))

# Create Seurat object with quality filtering
pbmc <- CreateSeuratObject(
  counts = raw_counts,
  project = PROJECT_NAME,
  min.cells = 3,      # Filter genes expressed in < 3 cells
  min.features = 200  # Filter cells with < 200 genes
)

message(paste("After initial filtering:", ncol(pbmc), "cells and", nrow(pbmc), "genes"))

# ==============================================================================
# 2. QUALITY CONTROL METRICS
# ==============================================================================

message("Calculating quality control metrics...")

# Calculate mitochondrial gene percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Store QC metrics before filtering
qc_metrics_before <- data.frame(
  n_cells = ncol(pbmc),
  n_genes = nrow(pbmc)
)

# Visualize QC metrics before filtering
qc_plot_before <- VlnPlot(
  pbmc,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.1
) + plot_annotation(title = "Quality Control Metrics - Before Filtering")

ggsave(
  file.path(FIGURES_DIR, "01_qc_metrics_before_filtering.png"),
  qc_plot_before,
  width = 15, height = 8, dpi = DPI, bg = "white"
)

# Create scatter plots for QC relationships
qc_scatter1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
qc_scatter2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

qc_scatter_combined <- qc_scatter1 + qc_scatter2
ggsave(
  file.path(FIGURES_DIR, "01_qc_scatter_plots.png"),
  qc_scatter_combined,
  width = 12, height = 6, dpi = DPI, bg = "white"
)

# ==============================================================================
# 3. CELL FILTERING
# ==============================================================================

message("Applying quality control filters...")

# Define filtering thresholds
MIN_GENES <- 200
MAX_GENES <- 5000
MAX_MT_PERCENT <- 15

# Apply filters
pbmc <- subset(
  pbmc,
  subset = nFeature_RNA > MIN_GENES & 
    nFeature_RNA < MAX_GENES & 
    percent.mt < MAX_MT_PERCENT
)

message(paste("After QC filtering:", ncol(pbmc), "cells and", nrow(pbmc), "genes"))

# Store QC metrics after filtering
qc_metrics_after <- data.frame(
  n_cells = ncol(pbmc),
  n_genes = nrow(pbmc),
  median_genes_per_cell = median(pbmc$nFeature_RNA),
  median_UMI_per_cell = median(pbmc$nCount_RNA),
  median_mt_percent = median(pbmc$percent.mt)
)

# Visualize QC metrics after filtering
qc_plot_after <- VlnPlot(
  pbmc,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.1
) + plot_annotation(title = "Quality Control Metrics - After Filtering")

ggsave(
  file.path(FIGURES_DIR, "01_qc_metrics_after_filtering.png"),
  qc_plot_after,
  width = 15, height = 8, dpi = DPI, bg = "white"
)

# ==============================================================================
# 4. NORMALIZATION
# ==============================================================================

message("Normalizing expression data...")

# Log-normalize the data
pbmc <- NormalizeData(
  pbmc,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# ==============================================================================
# 5. FEATURE SELECTION
# ==============================================================================

message("Identifying highly variable features...")

# Find highly variable features
pbmc <- FindVariableFeatures(
  pbmc,
  selection.method = "vst",
  nfeatures = 2000
)

# Get top 10 most variable genes
top10_variable <- head(VariableFeatures(pbmc), 10)
message(paste("Top 10 variable genes:", paste(top10_variable, collapse = ", ")))

# Plot variable features
variable_plot <- VariableFeaturePlot(pbmc)
variable_plot_labeled <- LabelPoints(
  plot = variable_plot,
  points = top10_variable,
  repel = TRUE,
  xnudge = 0,
  ynudge = 0
)

ggsave(
  file.path(FIGURES_DIR, "02_highly_variable_features.png"),
  variable_plot_labeled,
  width = 12, height = 8, dpi = DPI, bg = "white"
)

# ==============================================================================
# 6. DATA SCALING
# ==============================================================================

message("Scaling data...")

# Scale all genes
all_genes <- rownames(pbmc)
pbmc <- ScaleData(
  pbmc,
  features = all_genes
)

# ==============================================================================
# 7. PRINCIPAL COMPONENT ANALYSIS
# ==============================================================================

message("Running Principal Component Analysis...")

# Perform PCA
pbmc <- RunPCA(
  pbmc,
  features = VariableFeatures(object = pbmc),
  seed.use = RANDOM_SEED
)

# Visualize PCA results
pca_plot <- DimPlot(
  pbmc,
  reduction = "pca",
  dims = c(1, 2)
) + 
  ggtitle("PCA Plot") +
  theme(legend.position = "none")

ggsave(
  file.path(FIGURES_DIR, "03_pca_plot.png"),
  pca_plot,
  width = 8, height = 6, dpi = DPI, bg = "white"
)

# Create PCA heatmap
pca_heatmap <- DimHeatmap(
  pbmc,
  dims = 1:9,
  cells = 500,
  balanced = TRUE,
  fast = FALSE
)

ggsave(
  file.path(FIGURES_DIR, "03_pca_heatmap.png"),
  pca_heatmap,
  width = 12, height = 15, dpi = DPI, bg = "white"
)

# Create elbow plot for PC selection
elbow_plot <- ElbowPlot(pbmc, ndims = 20) +
  ggtitle("PCA Elbow Plot") +
  geom_vline(xintercept = 12, linetype = "dashed", color = "red") +
  annotate("text", x = 12, y = max(pbmc@reductions$pca@stdev[1:20]), 
           label = "Selected: PC12", hjust = -0.1, color = "red")

ggsave(
  file.path(FIGURES_DIR, "03_elbow_plot.png"),
  elbow_plot,
  width = 10, height = 6, dpi = DPI, bg = "white"
)

# ==============================================================================
# 8. CLUSTERING
# ==============================================================================

message("Performing clustering analysis...")

# Selected number of PCs based on elbow plot
N_PCS <- 12
RESOLUTION <- 0.5

# Find neighbors
pbmc <- FindNeighbors(
  pbmc,
  dims = 1:N_PCS,
  k.param = 20
)

# Find clusters
pbmc <- FindClusters(
  pbmc,
  resolution = RESOLUTION,
  random.seed = RANDOM_SEED
)

message(paste("Identified", length(unique(pbmc$seurat_clusters)), "clusters"))

# ==============================================================================
# 9. DIMENSIONALITY REDUCTION (UMAP)
# ==============================================================================

message("Running UMAP for visualization...")

# Run UMAP
pbmc <- RunUMAP(
  pbmc,
  dims = 1:N_PCS,
  seed.use = RANDOM_SEED
)

# Create UMAP plot
umap_plot <- DimPlot(
  pbmc,
  reduction = "umap",
  label = TRUE,
  label.size = 6,
  repel = TRUE
) + 
  ggtitle("UMAP Clustering Results") +
  theme(legend.position = "right")

ggsave(
  file.path(FIGURES_DIR, "04_umap_clusters.png"),
  umap_plot,
  width = 12, height = 8, dpi = DPI, bg = "white"
)

# Create UMAP plot without labels for cleaner look
umap_plot_clean <- DimPlot(
  pbmc,
  reduction = "umap",
  label = FALSE
) + 
  ggtitle("UMAP Clustering Results") +
  theme(legend.position = "right")

ggsave(
  file.path(FIGURES_DIR, "04_umap_clusters_clean.png"),
  umap_plot_clean,
  width = 12, height = 8, dpi = DPI, bg = "white"
)

# ==============================================================================
# 10. SAVE PROCESSED DATA
# ==============================================================================

message("Saving processed Seurat object...")

# Save the processed Seurat object
saveRDS(pbmc, file = file.path(RESULTS_DIR, "processed_pbmc.rds"))

# ==============================================================================
# 11. SUMMARY REPORT
# ==============================================================================

message("Generating summary report...")

# Create summary statistics
summary_stats <- data.frame(
  Metric = c(
    "Initial cells",
    "Initial genes",
    "Cells after QC",
    "Genes after QC",
    "Highly variable features",
    "Principal components used",
    "Clusters identified",
    "Clustering resolution"
  ),
  Value = c(
    qc_metrics_before$n_cells,
    qc_metrics_before$n_genes,
    qc_metrics_after$n_cells,
    qc_metrics_after$n_genes,
    length(VariableFeatures(pbmc)),
    N_PCS,
    length(unique(pbmc$seurat_clusters)),
    RESOLUTION
  )
)

# Print summary
message("\n=== PREPROCESSING SUMMARY ===")
for (i in 1:nrow(summary_stats)) {
  message(paste(summary_stats$Metric[i], ":", summary_stats$Value[i]))
}

# Save summary to file
write.csv(summary_stats, file = file.path(RESULTS_DIR, "preprocessing_summary.csv"), row.names = FALSE)

# Save QC metrics comparison
qc_comparison <- data.frame(
  Stage = c("Before QC", "After QC"),
  Cells = c(qc_metrics_before$n_cells, qc_metrics_after$n_cells),
  Genes = c(qc_metrics_before$n_genes, qc_metrics_after$n_genes),
  Median_Genes_per_Cell = c(qc_metrics_before$median_genes_per_cell, qc_metrics_after$median_genes_per_cell),
  Median_UMI_per_Cell = c(qc_metrics_before$median_UMI_per_cell, qc_metrics_after$median_UMI_per_cell),
  Median_MT_Percent = c(qc_metrics_before$median_mt_percent, qc_metrics_after$median_mt_percent)
)

write.csv(qc_comparison, file = file.path(RESULTS_DIR, "qc_metrics_comparison.csv"), row.names = FALSE)

# Save cluster information
cluster_info <- data.frame(
  Cluster = names(table(pbmc$seurat_clusters)),
  Cell_Count = as.numeric(table(pbmc$seurat_clusters))
)

write.csv(cluster_info, file = file.path(RESULTS_DIR, "cluster_summary.csv"), row.names = FALSE)

# Save highly variable genes
variable_genes_df <- data.frame(
  Gene = VariableFeatures(pbmc)[1:50],  # Top 50 variable genes
  Rank = 1:50
)

write.csv(variable_genes_df, file = file.path(RESULTS_DIR, "top_variable_genes.csv"), row.names = FALSE)

message("\nPreprocessing pipeline completed successfully!")
message(paste("Results saved to:", getwd()))
message("Generated files:")
message(paste("- Processed Seurat object:", file.path(RESULTS_DIR, "processed_pbmcect.rds")))
message(paste("- Summary statistics:", file.path(RESULTS_DIR, "preprocessing_summary.csv")))
message(paste("- QC comparison:", file.path(RESULTS_DIR, "qc_metrics_comparison.csv")))
message(paste("- Cluster information:", file.path(RESULTS_DIR, "cluster_summary.csv")))
message(paste("- Top variable genes:", file.path(RESULTS_DIR, "top_variable_genes.csv")))
message(paste("- Figures in", FIGURES_DIR, "directory"))

# ==============================================================================
# END OF SCRIPT
# ==============================================================================
