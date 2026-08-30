# ==============================================================================
# Figure 4A-D: Trajectory, pseudotime and branch-specific regulatory analysis
#
# Description:
# Generates Monocle trajectory and pseudotime visualizations for integrated
# EC and beta-cell populations, branch-dependent gene-expression profiles,
# and differential motif activity between EC and beta-cell branches.
#
# Inputs:
#   Monocle trajectory object containing integrated EC and beta-cell populations
#   Integrated EC/beta-cell multiome dataset with Monocle trajectory-state annotations
#
# Outputs:
#   figures/Figure4/
#
# Environment:
#   R package versions are recorded in renv.lock.
# ==============================================================================

# ==============================================================================
# Packages and settings
# ==============================================================================

library(monocle)
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnhancedVolcano)
set.seed(1234)

# ==============================================================================
# Input data and output directory
# ==============================================================================

HSMM_myo <- readRDS(
  "data/processed/ec_beta_hsmm_monocle_obj.rds"
)

integrated_for_monocle <- readRDS(
  "data/processed/integrated_for_monocle.rds"
)

output_dir <- "figures/Figure4"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Shared plotting settings
# ==============================================================================

trajectory_colors <- c(
  "Beta cells-1" = "#FF61CC",
  "Hi INS Beta cells-2" = "#C77CFF",
  "EC cells-1" = "#F8766D",
  "EC cells-4" = "#E68613",
  "EC cells-5" = "#ABA300",
  "EC cells-2" = "#00BCD7",
  "EC cells-3" = "#7CAE00"
)

# ==============================================================================
# Figure 4A: Monocle trajectory and pseudotime
# ==============================================================================

trj_cluster<- plot_cell_trajectory(HSMM_myo, color_by = "celltype") & scale_color_manual(values = trajectory_colors)

trj_pseudotime<- plot_cell_trajectory(HSMM_myo, color_by = "Pseudotime" )+theme(legend.position =  c(0.9, 0.22),legend.text = element_text(size=12),legend.title = element_text(size=12),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=1),axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),axis.ticks = element_line(size = 1),axis.text.x=element_text(size = 14,colour = 'black'),axis.text.y=element_text(size = 14,colour = 'black'))
ggsave(
  plot = trj_cluster,
  filename = file.path(output_dir, "Figure4A_trajectory_celltype.pdf"),
  width = 6,
  height = 6
)

ggsave(
  plot = trj_pseudotime,
  filename = file.path(output_dir, "Figure4A_trajectory_pseudotime.pdf"),
  width = 6,
  height = 6
)

# ==============================================================================
# Figure 4B: Branch-dependent gene-expression heatmap
# ==============================================================================
BEAM_res <- BEAM(HSMM_myo, branch_point = 2, cores = 4)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
BEAM_res %>% arrange(pval) %>% head(100) -> beam_top100
heatmap_branch<- plot_genes_branched_heatmap(HSMM_myo[row.names(beam_top100),],
                                          branch_point = 2,
                                          num_clusters = 3,
                                          cores = 4,
                                          use_gene_short_name = T,
                                          show_rownames = T, return_heatmap=TRUE)
 ggsave(
  plot = heatmap_branch,
  filename = file.path(output_dir, "Figure4B_heatmap_monoclebranch.pdf"),
  width = 6,
  height = 6
)                                         
                                          
# ==============================================================================
# Figure 4C: Gene expression along pseudotime
# ==============================================================================

to_be_tested <- row.names(subset(fData(HSMM_myo),
gene_short_name %in% c( "EBF1","SLC30A8", "INS")))
cds_subset <-  HSMM_myo[to_be_tested,]

pseudotime_expression<- plot_genes_in_pseudotime(cds_subset, color_by = "celltype") & scale_color_manual(values = trajectory_colors)

test_genes <- row.names(subset(fData(HSMM_myo),
          gene_short_name %in% c("LMX1A", "LMX1B", "DDC")))

ec_branch_expression<- plot_genes_branched_pseudotime(HSMM_myo[test_genes,],
                       branch_point = 1,
                       color_by = "celltype",
                       ncol = 1) &scale_color_manual(values = trajectory_colors)
ggsave(
  plot = pseudotime_expression,
  filename = file.path(output_dir, "Figure4C_pseudotime_expression.pdf"),
  width = 6,
  height = 6
)

ggsave(
  plot = ec_branch_expression,
  filename = file.path(output_dir, "Figure4C_EC_branch_markers.pdf"),
  width = 6,
  height = 6
)
# ==============================================================================
# Figure 4D: Differential motif activity between trajectory branches
# ==============================================================================

DefaultAssay(integrated_for_monocle)<- 'ATAC'
integrated_for_monocle <- RunChromVAR(
  object = integrated_for_monocle,
  genome = BSgenome.Hsapiens.UCSC.hg38
)
DefaultAssay(integrated_for_monocle) <- 'chromvar'
differential.activity_monocle_branched <- FindMarkers(
  object = integrated_for_monocle,
  ident.1 = '3', #beta-cell branch
  ident.2 = '4', #ec cell branch
  only.pos = FALSE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)
DefaultAssay(integrated_for_monocle) <- 'ATAC'

differential.activity_monocle_branched$gene <- ConvertMotifID(integrated_for_monocle, id = rownames(differential.activity_monocle_branched))
differential.activity_monocle_branched<- differential.activity_monocle_branched %>% arrange(desc(avg_diff))
selectLab = c('MAF','MAFA', 'MAFG', 'NRL','MAFF','NEUROG2','NEUROD1','HOXA4', 'BACH1','HAND2')

branch_motif_volcano<- EnhancedVolcano(differential.activity_monocle_branched,
    lab = differential.activity_monocle_branched$gene, FCcutoff = 0.5,
    x = 'avg_diff', selectLab= selectLab,legendPosition = 'bottom',
    y = 'p_val',pointSize = 6.0,title='Motif Chromatin Accessibility',
                labSize = 6.0, drawConnectors=TRUE)
pdf(file.path(
    output_dir,
    "volcano_chromvar_motif_activ_beta_branch_vs_ec_top10_TFexpression.pdf"), width=10, height=10)
print(branch_motif_volcano)
dev.off()


# ==============================================================================
# Supplementary Figure 4: Motif enrichment in branch-specific accessible regions
#
# Description:
# Identifies differentially accessible regions between the beta-cell and EC
# trajectory branches and tests branch-specific accessible regions for
# transcription factor motif enrichment.
#
# Inputs:
#   Integrated EC/beta-cell multiome dataset with Monocle trajectory-state
#   annotations
#
# Outputs:
#   figures/SupplementaryFigure4/
#
# Environment:
#   R package versions are recorded in renv.lock.
# ==============================================================================

library(Seurat)
library(Signac)
library(ggplot2)
library(ggrepel)

set.seed(1234)

integrated_data <- readRDS(
  "data/processed/ec_beta_integrated_multiome.rds"
)

output_dir <- "figures/SupplementaryFigure4"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

DefaultAssay(integrated_data) <- "ATAC"

# Monocle trajectory states:
# state 3 = beta-cell branch
# state 4 = EC branch

# ==============================================================================
# Supplementary Figure 4A: Beta-branch motif enrichment
# ==============================================================================

beta_branch_da <- FindMarkers(
  object = integrated_data,
  ident.1 = "3",
  ident.2 = "4",
  only.pos = TRUE,
  test.use = "LR",
  min.pct = 0.05,
  latent.vars = "nCount_ATAC"
)

beta_branch_peaks <- rownames(
  beta_branch_da[beta_branch_da$p_val_adj < 0.05, ]
)

beta_branch_motifs <- FindMotifs(
  integrated_data,
  features = beta_branch_peaks
)

beta_branch_motifs <- beta_branch_motifs[
  beta_branch_motifs$p.adjust < 0.05,
]

beta_motif_plot <- ggplot(
  beta_branch_motifs,
  aes(x = -log10(pvalue), y = fold.enrichment)
) +
  geom_point(size = 5) +
  geom_text_repel(
    aes(label = motif.name),
    size = 5
  ) +
  theme_classic() +
  labs(
    title = "Motifs in DARs of the beta-cell branch",
    x = "-log10(p-value)",
    y = "Motif fold enrichment"
  )

ggsave(
  plot = beta_motif_plot,
  filename = file.path(
    output_dir,
    "beta_branch_DAR_motif_enrichment.pdf"
  ),
  width = 10,
  height = 10
)

# ==============================================================================
# Supplementary Figure 4B: EC-branch motif enrichment
# ==============================================================================

ec_branch_da <- FindMarkers(
  object = integrated_data,
  ident.1 = "4",
  ident.2 = "3",
  only.pos = TRUE,
  test.use = "LR",
  min.pct = 0.05,
  latent.vars = "nCount_ATAC"
)

ec_branch_peaks <- rownames(
  ec_branch_da[ec_branch_da$p_val_adj < 0.05, ]
)

ec_branch_motifs <- FindMotifs(
  integrated_data,
  features = ec_branch_peaks
)

ec_branch_motifs <- ec_branch_motifs[
  ec_branch_motifs$p.adjust < 0.05,
]

ec_motif_plot <- ggplot(
  ec_branch_motifs,
  aes(x = -log10(pvalue), y = fold.enrichment)
) +
  geom_point(size = 5) +
  geom_text_repel(
    aes(label = motif.name),
    size = 5
  ) +
  theme_classic() +
  labs(
    title = "Motifs in DARs of the EC branch",
    x = "-log10(p-value)",
    y = "Motif fold enrichment"
  )

ggsave(
  plot = ec_motif_plot,
  filename = file.path(
    output_dir,
    "EC_branch_DAR_motif_enrichment.pdf"
  ),
  width = 10,
  height = 10
)