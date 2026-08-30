# ==============================================================================
# Day 30 and Day 50 single-cell multiome processing and integration
#
# Description:
# Performs quality control, RNA and ATAC preprocessing, multimodal integration,
# weighted-nearest-neighbor analysis, clustering, and marker-based cell-type
# annotation of Day 30 and Day 50 single-cell multiome datasets.
#
# Inputs:
# data/raw/D30_filtered_feature_bc_matrix.h5
# data/raw/D30_atac_fragments.tsv.gz
# data/raw/D50_filtered_feature_bc_matrix.h5
# data/raw/D50_atac_fragments.tsv.gz
#
# Outputs:
# data/processed/integrated_D30_D50_multiome.rds
#
# Environment:
# R package versions are recorded in renv.lock.
# ==============================================================================


# ==============================================================================
# Packages and settings
# ==============================================================================

library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(ggplot2)

set.seed(1234)


# ==============================================================================
# Genome annotation
# ==============================================================================

annotation <- GetGRangesFromEnsDb(
  ensdb = EnsDb.Hsapiens.v86
)

genome(annotation) <- "hg38"
seqlevelsStyle(annotation) <- "UCSC"


# ==============================================================================
# Load Day 30 multiome data
# ==============================================================================

counts30 <- Read10X_h5(
  "data/raw/D30_filtered_feature_bc_matrix.h5"
)

fragpath30 <- "data/raw/D30_atac_fragments.tsv.gz"

d30_data <- CreateSeuratObject(
  counts = counts30$`Gene Expression`,
  assay = "RNA"
)

d30_data[["ATAC"]] <- CreateChromatinAssay(
  counts = counts30$Peaks,
  sep = c(":", "-"),
  fragments = fragpath30,
  annotation = annotation
)

d30_data$orig.ident <- "D30"

DefaultAssay(d30_data) <- "ATAC"
d30_data <- NucleosomeSignal(d30_data)
d30_data <- TSSEnrichment(d30_data)

DefaultAssay(d30_data) <- "RNA"
d30_data[["percent.mt"]] <- PercentageFeatureSet(
  d30_data,
  pattern = "^MT-"
)

d30_data$mitoRatio <- d30_data$percent.mt / 100


# ==============================================================================
# Load Day 50 multiome data
# ==============================================================================

counts50 <- Read10X_h5(
  "data/raw/D50_filtered_feature_bc_matrix.h5"
)

fragpath50 <- "data/raw/D50_atac_fragments.tsv.gz"

d50_data <- CreateSeuratObject(
  counts = counts50$`Gene Expression`,
  assay = "RNA"
)

d50_data[["ATAC"]] <- CreateChromatinAssay(
  counts = counts50$Peaks,
  sep = c(":", "-"),
  fragments = fragpath50,
  annotation = annotation
)

d50_data$orig.ident <- "D50"

DefaultAssay(d50_data) <- "ATAC"
d50_data <- NucleosomeSignal(d50_data)
d50_data <- TSSEnrichment(d50_data)

DefaultAssay(d50_data) <- "RNA"
d50_data[["percent.mt"]] <- PercentageFeatureSet(
  d50_data,
  pattern = "^MT-"
)

d50_data$mitoRatio <- d50_data$percent.mt / 100


# ==============================================================================
# Merge datasets and calculate QC metrics
# ==============================================================================

merged_seurat <- merge(
  x = d30_data,
  y = d50_data,
  add.cell.ids = c("D30", "D50")
)

merged_seurat[["percent.mt"]] <- PercentageFeatureSet(
  merged_seurat,
  pattern = "^MT-"
)

merged_seurat$mitoRatio <- merged_seurat$percent.mt / 100

merged_seurat$log10GenesPerUMI <-
  log10(merged_seurat$nFeature_RNA) /
  log10(merged_seurat$nCount_RNA)


# ==============================================================================
# Quality-control filtering
# ==============================================================================

data_filtered <- subset(
  merged_seurat,
  subset =
    nCount_ATAC >= 300 &
    nCount_RNA >= 300 &
    nFeature_RNA >= 100 &
    nFeature_ATAC >= 100 &
    log10GenesPerUMI >= 0.75 &
    percent.mt <= 50 &
    nucleosome_signal < 1 &
    TSS.enrichment > 1
)

# Cells retained after filtering:
# Day 30: 3,320
# Day 50: 3,902


# ==============================================================================
# Separate Day 30 and Day 50 datasets for integration
# ==============================================================================

d30_filtered <- subset(
  data_filtered,
  subset = orig.ident == "D30"
)

d50_filtered <- subset(
  data_filtered,
  subset = orig.ident == "D50"
)


# ==============================================================================
# RNA preprocessing
# ==============================================================================

d30_filtered <- SCTransform(
  d30_filtered,
  assay = "RNA"
)

d30_filtered <- FindVariableFeatures(d30_filtered)
d30_filtered <- RunPCA(d30_filtered)

d50_filtered <- SCTransform(
  d50_filtered,
  assay = "RNA"
)

d50_filtered <- FindVariableFeatures(d50_filtered)
d50_filtered <- RunPCA(d50_filtered)


# ==============================================================================
# RNA integration
# ==============================================================================

integration_features <- SelectIntegrationFeatures(
  object.list = list(d30_filtered, d50_filtered)
)

integration_objects <- PrepSCTIntegration(
  object.list = list(d30_filtered, d50_filtered),
  anchor.features = integration_features
)

anchors <- FindIntegrationAnchors(
  object.list = integration_objects,
  reduction = "rpca",
  dims = 1:20,
  normalization.method = "SCT",
  anchor.features = integration_features,
  assay = c("SCT", "SCT")
)

integratedRNA <- IntegrateData(
  anchorset = anchors,
  normalization.method = "SCT",
  dims = 1:20,
  new.assay.name = "integratedRNA"
)


# ==============================================================================
# ATAC dimensional reduction and integration
# ==============================================================================

integrated <- integratedRNA

DefaultAssay(integrated) <- "ATAC"

integrated <- FindTopFeatures(
  integrated,
  min.cutoff = 5
)

integrated <- RunTFIDF(integrated)
integrated <- RunSVD(integrated)

integrated <- IntegrateEmbeddings(
  anchorset = anchors,
  new.reduction.name = "integratedLSI",
  reductions = integrated[["lsi"]]
)


# ==============================================================================
# Weighted-nearest-neighbor analysis
# ==============================================================================

DefaultAssay(integrated) <- "integratedRNA"

integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated, verbose = FALSE)

integrated <- FindMultiModalNeighbors(
  object = integrated,
  reduction.list = list("pca", "integratedLSI"),
  dims.list = list(1:24, 2:24),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

integrated <- RunUMAP(
  object = integrated,
  nn.name = "weighted.nn",
  assay = "SCT",
  verbose = TRUE
)

integrated <- FindClusters(
  integrated,
  resolution = 0.06,
  graph.name = "wsnn",
  algorithm = 3,
  verbose = FALSE
)

DimPlot(
  integrated,
  reduction = "umap",
  label = TRUE,
  repel = TRUE
)


# ==============================================================================
# Cell-type annotation
# ==============================================================================

# Cell clusters were manually annotated based on expression of established
# pancreatic endocrine and enterochromaffin cell-type marker genes.
#
# Numerical cluster identities are not hard-coded here because cluster numbers
# can vary between clustering runs. Final cell-type annotations used for
# downstream analyses were assigned based on marker-expression profiles and
# stored in the metadata of the final processed Seurat object.

annotation_markers <- c(
  # Beta-cell markers
  "INS", "IAPP", "PDX1", "SLC30A8", "ISL1",

  # Alpha-cell markers
  "GCG", "ARX", "IRX2", "DPP4",

  # SST-cell marker
  "SST",

  # Enterochromaffin-cell markers
  "SLC18A1", "LMX1A", "LMX1B", "DDC"
)

DotPlot(
  integrated,
  features = annotation_markers
) +
  RotatedAxis()


# ==============================================================================
# Save processed integrated dataset
# ==============================================================================

saveRDS(
  integrated,
  "data/processed/integrated_D30_D50_multiome.rds"
)