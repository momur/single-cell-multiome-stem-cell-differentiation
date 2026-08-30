# ==============================================================================
# Monocle trajectory analysis: EC and beta-cell populations
#
# Description:
# Constructs a Monocle 2 CellDataSet from the integrated Day 30/Day 50
# Seurat object and performs trajectory inference and pseudotime ordering
# for EC and beta-cell populations.
#
# Input:
# data/processed/integrated_data_input_formonocle.rds
#
# Output:
# data/processed/ec_beta_hsmm_monocle_obj.rds
#
# Environment:
# R package versions are recorded in renv.lock.
# ==============================================================================


# ==============================================================================
# Packages and settings
# ==============================================================================

library(Seurat)
library(monocle)
library(dplyr)

set.seed(1234)


# ==============================================================================
# Input data
# ==============================================================================

integrated <- readRDS(
  "data/processed/integrated_data_input_formonocle.rds"
)


# ==============================================================================
# Subset EC and beta-cell populations
# ==============================================================================

ec_beta_cells <- subset(
  integrated,
  idents = c(
    "EC cells-1",
    "EC cells-2",
    "EC cells-3",
    "EC cells-4",
    "EC cells-5",
    "Beta cells-1",
    "Hi INS Beta cells-2"
  )
)


# ==============================================================================
# Construct Monocle CellDataSet
# ==============================================================================

expression_matrix <- as(
  as.matrix(ec_beta_cells@assays$RNA@counts),
  "sparseMatrix"
)

phenotype_data <- new(
  "AnnotatedDataFrame",
  data = ec_beta_cells@meta.data
)

feature_data <- data.frame(
  gene_short_name = rownames(expression_matrix),
  row.names = rownames(expression_matrix)
)

feature_data <- new(
  "AnnotatedDataFrame",
  data = feature_data
)

HSMM_myo <- newCellDataSet(
  expression_matrix,
  phenoData = phenotype_data,
  featureData = feature_data,
  lowerDetectionLimit = 0.5,
  expressionFamily = negbinomial.size()
)


# ==============================================================================
# Estimate expression parameters and detect expressed genes
# ==============================================================================

HSMM_myo <- estimateSizeFactors(HSMM_myo)

HSMM_myo <- estimateDispersions(
  HSMM_myo,
  cores = 4
)

HSMM_myo <- detectGenes(
  HSMM_myo,
  min_expr = 0.1
)

expressed_genes <- rownames(
  subset(
    fData(HSMM_myo),
    num_cells_expressed >= 50
  )
)


# ==============================================================================
# Identify ordering genes
# ==============================================================================

clustering_DEG_genes <- differentialGeneTest(
  HSMM_myo[expressed_genes, ],
  fullModelFormulaStr = "~celltype",
  cores = 3
)

HSMM_ordering_genes <- rownames(
  clustering_DEG_genes[
    order(clustering_DEG_genes$qval),
  ]
)[1:1000]

HSMM_myo <- setOrderingFilter(
  HSMM_myo,
  ordering_genes = HSMM_ordering_genes
)


# ==============================================================================
# DDRTree trajectory inference
# ==============================================================================

HSMM_myo <- reduceDimension(
  HSMM_myo,
  method = "DDRTree"
)

HSMM_myo <- orderCells(HSMM_myo)


# ==============================================================================
# Set EC cells-5-containing state as trajectory root
# ==============================================================================

GM_state <- function(cds) {

  if (length(unique(pData(cds)$State)) > 1) {

    T0_counts <- table(
      pData(cds)$State,
      pData(cds)$celltype
    )[,"EC cells-5"]

    return(
      as.numeric(
        names(T0_counts)[
          which(T0_counts == max(T0_counts))
        ]
      )
    )

  } else {

    return(1)

  }
}

HSMM_myo <- orderCells(
  HSMM_myo,
  root_state = GM_state(HSMM_myo)
)


# ==============================================================================
# Save final trajectory object
# ==============================================================================

saveRDS(
  HSMM_myo,
  "data/processed/ec_beta_hsmm_monocle_obj.rds"
)