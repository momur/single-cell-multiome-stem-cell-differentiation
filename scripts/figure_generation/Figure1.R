# ==============================================================================
# Figure 1: Day 30 and Day 50 single-cell multiome characterization
#
# Description:
# Generates RNA, ATAC and WNN UMAPs, marker-expression plots and ATAC
# coverage plots for Day 30 and Day 50 samples.
#
# Inputs:
#   Processed and cell-type-annotated Day 30 multiome dataset
#   Processed and cell-type-annotated Day 50 multiome dataset
#
# Outputs:
#   figures/Figure1/
#
# Environment:
#   R package versions are recorded in renv.lock.
# ==============================================================================

# ==============================================================================
# Packages and settings
# ==============================================================================

library(Signac)
library(Seurat)
library(ggplot2)

set.seed(1234)
# ==============================================================================
# Input data and output directory
# ==============================================================================

my_cols <- c('Beta cells'='#FF61CC',
             'Alpha-Beta cells'='#00A9FF',
             'Polyhormonal cells'='#8494FF',
             'EC cells'='#E68613',
             'SST cells'='#0CB702',
             'Alpha cells'='#00C19A',
             'UNK'= '#d3d3d3')

d30_input <- readRDS(
  "data/processed/day30_annotated_multiome.rds"
)

d50_input <- readRDS(
  "data/processed/day50_annotated_multiome.rds"
)

output_dir <- "figures/Figure1"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Day 30
# ==============================================================================                        
d30_input_data2 <- RenameIdents(
  object = d30_input,
  "Beta cells-1" = "Beta cells",
  "Hi INS Beta cells-2" = "Beta cells",
  "Alpha-Beta cells" = "Alpha-Beta cells",
  "Alpha cells" = "Alpha cells",
  "SST+ cells" = "SST cells",
  "Polyhormonal cells" = "Polyhormonal cells",
  "EC cells-1" = "EC cells",
  "EC cells-2" = "EC cells",
  "EC cells-3" = "EC cells",
  "EC cells-4" = "EC cells",
  "UNK" = "UNK"
)
# Figure 1B: Day 30 RNA, ATAC and WNN UMAPs
d30_umap_rna <- DimPlot(d30_input_data2, reduction = "umap.rna", label=T,label.box =T,repel = T,pt.size = 1, cols = my_cols) + ggtitle('RNA only') + NoLegend() +
theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=1))
d30_umap_rna
ggsave(plot = d30_umap_rna, filename = file.path(output_dir,"D30_umap_w_label_boxes_RNA_four_clusters.pdf"), width = 6, height = 6)

d30_umap_atac <- DimPlot(d30_input_data2, reduction = "umap.atac", label=T,label.box =T,repel = T,pt.size = 1, cols = my_cols) + ggtitle('ATAC only') + NoLegend()+
theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=1))
d30_umap_atac
ggsave(plot = d30_umap_atac, filename = file.path(output_dir, "D30_umap_w_label_boxes_atac_title_four_clusters.pdf"), width = 6, height = 6)

d30_umap_wnn <- DimPlot(d30_input_data2, reduction = "wnn.umap", label=T,label.box =T,repel = T,pt.size = 1, cols = my_cols) + ggtitle('Day 30-Combined modalities')+ theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=1))
d30_umap_wnn
ggsave(plot = d30_umap_wnn, filename = file.path(output_dir,"D30_umap_w_label_boxes_combined_title_four_clusters.pdf"), width = 6, height = 6)

# Figure 1C: Day 30 marker expression
d30_dotplot<- DotPlot(d30_input_data2, features = c("INS", "eGFP", "SST", "LMX1A", "GCG"),  assay = 'SCT')
pdf(file.path(output_dir,'day30_hormoneexpression_dotplot.pdf'), width=6, height=6)
print(d30_dotplot)
dev.off()

# optional alternative for Figure 1C - not included in the current version of the manuscript figures 
day30_beta_violin<- VlnPlot(d30_input_data2, features = c("INS", "SLC30A8", "ERO1B", "NKX2-2"), ncol = 4, pt.size=0,assay = 'SCT',cols = my_cols)
ggsave(day30_beta_violin, filename = file.path(output_dir,"Violinplot3_hormones_maturation_mark_Day30_four_clusters_beta_markers.jpg"), width = 16, height = 4)

day30_alpha_violin<- VlnPlot(d30_input_data2, features = c("GCG", "GC", "LDB2", "ARX"), ncol = 4, pt.size=0,assay = 'SCT',cols = my_cols)
ggsave(day30_alpha_violin, filename = file.path(output_dir,"Violinplot3_hormones_maturation_mark_Day30_four_clusters_alpha.jpg"), width = 16, height = 4)

day30_ec_violin<- VlnPlot(d30_input_data2, features = c("LMX1A","DDC","TPH1", "MNX1"), ncol = 4, pt.size=0,assay = 'SCT',cols = my_cols)
ggsave(day30_ec_violin, filename = file.path(output_dir,"Violinplot3_hormones_maturation_mark_Day30_four_clusters_EC.jpg"), width = 16, height = 4)


# Figure 1D: ATAC coverage of hormones on selected cell types
DefaultAssay(d30_input_data2) <- "ATAC"

cell_list_d30 <- WhichCells(
  d30_input_data2,
  idents = c(
    "Beta cells",
    "Alpha cells",
    "SST cells",
    "EC cells"
  )
)

d30_four_celltypes <- d30_input_data2[, cell_list_d30]

d30_coverage<- CoveragePlot(object = d30_four_celltypes,  region = c("INS", "GCG", "SST"), extend.upstream = 500,  extend.downstream = 500)+ theme(text =element_text(size=12),axis.title.x = element_text(size=12), axis.title.x.top = element_text(size=12),
axis.text=element_text(size=12),axis.title=element_text(size=12,face="bold"),axis.ticks = element_line(size = 1),axis.text.x=element_text(size = 12,colour = 'black'),axis.text.y=element_text(size = 12,colour = 'black')) & scale_fill_manual(values = c("#FF61CC", "#00C19A", "#0CB702", "#E68613"))
pdf(file.path(output_dir,"hormones_coverageplot_d30only_fourclusters_fig1manuscript.pdf"), height = 6, width = 8) 
print(d30_coverage)
dev.off()

# ==============================================================================
# Day 50
# ==============================================================================
d50_input_fourcluster <- RenameIdents(object=d50_input, 
                           "Beta cells" = "Beta cells",
                           "Hi INS Beta cells" = "Beta cells",
                           "Alpha-Beta cells"="Alpha-Beta cells",
                           "Alpha cells"= "Alpha cells",
                           "SST cells" = "SST cells", 
                           "Polyhormonal cells"="Polyhormonal cells",
                           "EC cells-1"="EC cells",
                           "EC cells-2"="EC cells",
                           "EC cells-3"="EC cells",
                           "EC cells-4"="EC cells",
                           "UNK"= "UNK")

# Figure 1E: day 50 RNA, ATAC and WNN UMAPs
d50_umap_rna <- DimPlot(d50_input_fourcluster, reduction = "umap.rna", label=T,label.box =T,repel = T,pt.size = 1, cols = my_cols) + ggtitle('RNA only') + NoLegend() +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=1))
d50_umap_rna
ggsave(plot = d50_umap_rna, filename= file.path(output_dir,"Day50_umap_w_label_boxes_RNA_four_clusters.pdf"), width = 6, height = 6)

d50_umap_atac <- DimPlot(d50_input_fourcluster, reduction = "umap.atac", label=T,label.box =T,repel = T,pt.size = 1, cols = my_cols) + ggtitle('ATAC only') + NoLegend()+
theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=1))
d50_umap_atac
ggsave(plot = d50_umap_atac, filename = file.path(output_dir,"Day50_umap_w_label_boxes_atac_title_four_clusters.pdf"), width = 6, height = 6)

d50_umap_wnn <- DimPlot(d50_input_fourcluster, reduction = "wnn.umap", label=T,label.box =T,repel = T,pt.size = 1, cols = my_cols) + ggtitle('Day 50-Combined modalities')+
theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=1))
d50_umap_wnn
ggsave(plot = d50_umap_wnn, filename = file.path(output_dir,"Day50_umap_w_label_boxes_combined_title_four_clusters.pdf"), width = 6, height = 6)

# Figure 1F: Day 50 hormone expression plots 
d50_dotplot<- DotPlot(d50_input_fourcluster, features = c("INS", "eGFP", "SST", "LMX1A", "GCG"),  assay = 'SCT')
pdf(file.path(output_dir,"day50_hormoneexpression_dotplot.pdf"), width=6, height=6)
print(d50_dotplot)
dev.off()

# optional alternative for Figure 1F - not included in the current version of the manuscript figures 
d50_beta_violin<- VlnPlot(d50_input_fourcluster, features = c("INS", "SLC30A8", "ERO1B", "NKX6-1"), ncol = 4, pt.size=0,assay = 'SCT',cols = my_cols)
ggsave(d50_beta_violin, filename = file.path(output_dir,"Violinplot3_hormones_maturation_mark_Day50_four_clusters_four_beta_markers_horizontal.jpg"), width = 16, height = 4)

d50_alpha_violin<- VlnPlot(d50_input_fourcluster, features = c("GCG", "LDB2", "IRX2", "DPP4"), ncol = 4, pt.size=0,assay = 'SCT',cols = my_cols)
ggsave(d50_alpha_violin, filename = file.path(output_dir,"Violinplot3_hormones_maturation_mark_Day50_four_clusters_four_alpha_markers_horizontal.jpg"), width = 16, height = 4)

d50_ec_violin<- VlnPlot(d50_input_fourcluster, features = c("LMX1A", "MNX1", "DDC","TPH1"), ncol = 4, pt.size=0,assay = 'SCT',cols = my_cols)
ggsave(d50_ec_violin, filename = file.path(output_dir, "Violinplot3_hormones_maturation_mark_Day50_four_clusters_four_EC_markers_horizontal.jpg"), width = 16, height = 4)

# Figure 1G: Day 50 ATAC coverage of hormones on selected cell types
DefaultAssay(d50_input_fourcluster) <- "ATAC"

cell_list_d50 <- WhichCells(
  d50_input_fourcluster,
  idents = c(
    "Beta cells",
    "Alpha cells",
    "SST cells",
    "EC cells"
  )
)

d50_four_celltypes <- d50_input_fourcluster[, cell_list_d50]
d50_coverage <- CoveragePlot(object = d50_four_celltypes,  region = c("INS", "GCG", "SST"), extend.upstream = 500,  extend.downstream = 500)+ theme(text =element_text(size=12),axis.title.x = element_text(size=12), axis.title.x.top = element_text(size=12),
axis.text=element_text(size=12),axis.title=element_text(size=12,face="bold"), axis.ticks = element_line(size = 1),axis.text.x=element_text(size = 12,colour = 'black'),axis.text.y=element_text(size = 12,colour = 'black')) & scale_fill_manual(values = c("#FF61CC", "#00C19A", "#0CB702", "#E68613"))
pdf(file.path(output_dir,"hormones_coverageplot_day50only_fourclusters_fig1manuscript.pdf"), height = 6, width = 8)
print(d50_coverage)
dev.off()







                           