# ==============================================================================
# Figure 4E: Shared transcription factors across cell states from CellOracle GRN analysis
#
# Description:
# Generates an UpSet plot showing overlap among the top cell-type-specific
# transcription factors across EC and beta-cell states.
#
# Inputs:
#   data/processed/TF_clusters.csv
#
# Outputs:
#   figures/Figure4/Figure4E_TF_overlap_upset.pdf
#
# Environment:
#   R package versions are recorded in renv.lock.
# ==============================================================================
# ==============================================================================
# Packages and settings
# ==============================================================================

library(dplyr)
library(readr)
library(UpSetR)


# ==============================================================================
# Input data and output directory
# ==============================================================================

tf_clusters <- read_csv(
  "data/processed/TF_clusters.csv"
)

output_dir <- "figures/Figure4"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# Figure 4E: UpSet plot of shared transcription factors
# ==============================================================================

# Create a named list of unique TFs for each cell state
tf_lists <- tf_clusters %>%
  split(.$cluster) %>%
  lapply(function(x) unique(x$TF))

# Define the order of cell states in the plot
celltype_order <- c(
  "EC_cells_5",
  "EC_cells_4",
  "EC_cells_1",
  "EC_cells_2",
  "EC_cells_3",
  "Beta_cells_1",
  "HiINS_Beta_cells_2"
)

# Generate and save UpSet plot
pdf(
  file.path(output_dir, "Figure4E_TF_overlap_upset.pdf"),
  width = 9,
  height = 6
)

upset(
  fromList(tf_lists),
  nsets = length(tf_lists),
  order.by = "freq",
  sets = celltype_order,
  text.scale = 1.2,
  point.size = 3,
  mb.ratio = c(0.2, 0.8)
)

dev.off()