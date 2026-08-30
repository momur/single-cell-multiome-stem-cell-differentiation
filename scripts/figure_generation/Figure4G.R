# ==============================================================================
# Figure 4G: Cell-type-specific regulon activity plot generation
#
# Description:
# Loads SCENIC regulon activity results from a loom file, calculates regulon
# specificity scores (RSS) across EC and beta-cell states, and generates the
# Figure 4G regulon specificity plot.
#
# The script also generates the supplementary heatmap of selected regulons.
#
# Inputs:
#   SCENIC AUCell loom file for integrated EC and beta-cell populations
#
# Outputs:
#   figures/Figure4/
#   figures/Supplementary/
#
# Environment:
#   R package versions are recorded in renv.lock.
# ==============================================================================


# ==============================================================================
# Packages and settings
# ==============================================================================

library(SCENIC)
library(AUCell)
library(loomR)
library(SCopeLoomR)

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

library(ComplexHeatmap)
library(circlize)
library(grid)


# ==============================================================================
# Input data and output directories
# ==============================================================================

scenic_loom_path <- "data/processed/ec_and_beta_integ_AUCell_singleTF_1.loom"

output_dir <- "figures/Figure4"
supp_output_dir <- "figures/Supplementary"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_output_dir, recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# Load SCENIC regulon activity
# ==============================================================================

loom <- open_loom(scenic_loom_path)

regulon_auc <- get_regulons_AUC(
  loom,
  column.attr.name = "RegulonsAUC"
)

cell_clusters <- get_cell_annotation(loom)

close_loom(loom)


# Remove duplicated extended regulons
regulon_auc <- regulon_auc[
  onlyNonDuplicatedExtended(rownames(regulon_auc)),
]


# ==============================================================================
# Calculate average regulon activity by cell type
# ==============================================================================

selected_resolution <- "celltype"

cells_per_cluster <- split(
  rownames(cell_clusters),
  cell_clusters[, selected_resolution]
)

regulon_activity_by_celltype <- sapply(
  cells_per_cluster,
  function(cells) {
    rowMeans(getAUC(regulon_auc)[, cells, drop = FALSE])
  }
)

# Remove cell-type columns containing only missing values
regulon_activity_by_celltype <- regulon_activity_by_celltype[
  ,
  colSums(is.na(regulon_activity_by_celltype)) <
    nrow(regulon_activity_by_celltype),
  drop = FALSE
]


# ==============================================================================
# Calculate regulon specificity scores
# ==============================================================================

rss <- calcRSS(
  AUC = getAUC(regulon_auc),
  cellAnnotation = cell_clusters[
    colnames(regulon_auc),
    selected_resolution
  ]
)


# ==============================================================================
# Figure 4G: Regulon specificity across EC and beta-cell states
# ==============================================================================

celltype_order <- c(
  "Hi INS Beta cells-2",
  "Beta cells-1",
  "EC cells-3",
  "EC cells-2",
  "EC cells-1",
  "EC cells-4",
  "EC cells-5"
)

cell_colors <- c(
  "Hi INS Beta cells-2" = "#C77CFF",
  "Beta cells-1" = "#FF61CC",
  "EC cells-3" = "#00BE67",
  "EC cells-2" = "#00BCD7",
  "EC cells-1" = "#F8766D",
  "EC cells-4" = "#E68613",
  "EC cells-5" = "#ABA300"
)

cell_colors_sub <- cell_colors[celltype_order]


tfs_to_plot <- c(
  "EBF1(+)",
  "PDX1(+)",
  "STAT1(+)",
  "TCF4(+)",
  "FOXA2(+)",
  "MNX1(+)",
  "ESR1(+)",
  "NR3C1(+)",
  "IRF1(+)",
  "MEF2C(+)",
  "PLAGL1(+)",
  "NR1H4(+)",
  "KLF5(+)",
  "PPARGC1A(+)",
  "ZBTB7C(+)",
  "TCF7L2(+)",
  "NPAS3(+)",
  "ETS1(+)"
)


# Convert RSS matrix to long format
rss_long <- as.data.frame(rss) %>%
  rownames_to_column("TF") %>%
  filter(TF %in% tfs_to_plot) %>%
  pivot_longer(
    -TF,
    names_to = "cell_type",
    values_to = "RSS"
  ) %>%
  filter(cell_type %in% celltype_order) %>%
  mutate(
    TF = factor(TF, levels = tfs_to_plot),
    cell_type = factor(
      cell_type,
      levels = celltype_order
    )
  )


regulon_specificity_plot <- ggplot(
  rss_long,
  aes(
    x = cell_type,
    y = RSS,
    color = cell_type
  )
) +
  geom_point(size = 2) +
  scale_color_manual(values = cell_colors_sub) +
  facet_wrap(
    ~ TF,
    scales = "free_y",
    ncol = 5
  ) +
  labs(
    x = "Cell state",
    y = "Regulon Specificity Score (RSS)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(
      face = "bold.italic",
      size = 8
    ),
    strip.background = element_rect(
      fill = "white",
      color = "black"
    ),
    legend.title = element_blank(),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)
  ) +
  guides(
    color = guide_legend(nrow = 1)
  )


ggsave(
  plot = regulon_specificity_plot,
  filename = file.path(
    output_dir,
    "Figure4G_regulon_specificity.pdf"
  ),
  width = 8,
  height = 4
)


# ==============================================================================
# Supplementary figure: Regulon activity heatmap
# ==============================================================================

# NOTE:
# `topRegulators` must be generated or loaded by the upstream SCENIC workflow
# before this section can be reproduced.


top25_regulons <- topRegulators %>%
  group_by(CellType) %>%
  arrange(
    desc(RelativeActivity),
    .by_group = TRUE
  ) %>%
  slice_head(n = 25) %>%
  ungroup() %>%
  pull(Regulon) %>%
  unique()


# Subset average regulon activity
regulon_matrix <- regulon_activity_by_celltype[
  rownames(regulon_activity_by_celltype) %in% top25_regulons,
  ,
  drop = FALSE
]


# Reorder cell types
heatmap_celltype_order <- intersect(
  celltype_order,
  colnames(regulon_matrix)
)

regulon_matrix <- regulon_matrix[
  ,
  heatmap_celltype_order,
  drop = FALSE
]


# Remove regulons with zero variance
regulon_matrix <- regulon_matrix[
  apply(
    regulon_matrix,
    1,
    sd,
    na.rm = TRUE
  ) > 0,
  ,
  drop = FALSE
]


# Scale activity by regulon
regulon_matrix_scaled <- t(
  scale(t(regulon_matrix))
)


# Remove non-finite values
regulon_matrix_scaled <- regulon_matrix_scaled[
  apply(
    regulon_matrix_scaled,
    1,
    function(x) all(is.finite(x))
  ),
  apply(
    regulon_matrix_scaled,
    2,
    function(x) all(is.finite(x))
  ),
  drop = FALSE
]


heatmap_colors <- circlize::colorRamp2(
  c(-2, 0, 2),
  hcl_palette = "Tropic"
)


regulon_heatmap <- ComplexHeatmap::Heatmap(
  regulon_matrix_scaled,
  cluster_columns = FALSE,
  column_order = heatmap_celltype_order,
  name = "Regulon activity",
  col = heatmap_colors,
  row_names_gp = grid::gpar(fontsize = 3),
  column_names_gp = grid::gpar(fontsize = 5),
  column_names_rot = 45
)


pdf(
  file.path(
    supp_output_dir,
    "top25_regulator_heatmap_ecandbeta.pdf"
  ),
  height = 6,
  width = 4
)

ComplexHeatmap::draw(regulon_heatmap)

dev.off()