# ==============================================================================
# Figure 3: Integrated Day 30 and Day 50 multiome analysis
#
# Description:
# Generates visualizations from the integrated Day 30 and Day 50 single-cell
# multiome dataset, including integrated UMAPs, EC-subcluster marker expression,
# transcription factor expression and motif activity, cell-type marker heatmaps,
# and beta-cell marker expression.
#
# Inputs:
#   Integrated and cell-type-annotated Day 30/Day 50 multiome dataset

# Outputs:
#   figures/Figure3/
#
# Environment:
#   R package versions are recorded in renv.lock.
# ==============================================================================

# ==============================================================================
# Packages and settings
# ==============================================================================

library(Seurat)
library(Signac)
library(ggplot2)
library(circlize)
library(ComplexHeatmap)
library(gridExtra)
set.seed(1234)

my_cols_int <- c('Beta cells'='#FF61CC',
             'Hi INS Beta cells'='#C77CFF',
             'Alpha-Beta cells'='#00A9FF',
             'Polyhormonal cells'='#8494FF',
             'EC cells-1'='#F8766D',
             'EC cells-4'='#E68613',
             'EC cells-5'='#ABA300',
             'EC cells-2'='#00BCD7',
             'EC cells-3'='#00BE67',
             'SST cells'='#0CB702',
             'Alpha cells'='#00C19A',
             'UNK'= '#d3d3d3')
# ==============================================================================
# Input data and output directory
# ==============================================================================

integ_data <- readRDS(
  "data/processed/day30_day50_integrated_annotated_multiome.rds"
)

output_dir <- "figures/Figure3"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

#Figure 3A
integ_data <- RenameIdents(object=integ_data, 
                           "Beta cells-1" = "Beta cells",
                           "Hi INS Beta cells-2" = "Hi INS Beta cells",
                           "Alpha-Beta cells 1"="Alpha-Beta cells",
                           "Alpha cells"= "Alpha cells",
                           "SST+ cells" = "SST cells", 
                           "Polyhormonal cells"="Polyhormonal cells",
                           "EC cells-1"="EC cells-1",
                           "EC cells-2"="EC cells-2",
                           "EC cells-3"="EC cells-3",
                           "EC cells-4"="EC cells-4",
                           "EC cells-5"="EC cells-5",
                           "UNK"= "UNK")
                      
# ==============================================================================
# Figure 3A: UMAP of integrated datasets of Day 30 and Day 50
# ==============================================================================

integrated_umap<- DimPlot(integ_data, reduction = "umap", label=T,label.box =T,repel = T,pt.size = 1, cols = my_cols_int )+
theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=1)) +ggtitle('Integrated Data Day 30 & Day 50')
ggsave(plot = integrated_umap, filename = file.path(output_dir,"Integrated_umap_w_label_boxes_legendposi_changed_fullframe.pdf"), width = 9, height = 8)       
integ_sep <- DimPlot(
  integ_data,
  reduction = "umap",
  repel = TRUE,
  pt.size = 0.25,
  split.by = "orig.ident"
)

integ_plot <- integ_sep +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour = "black", linewidth = 1),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 10, face = "bold"),
    axis.ticks = element_line(linewidth = 1),
    axis.text.x = element_text(size = 10, colour = "black"),
    axis.text.y = element_text(size = 10, colour = "black")
  )
ggsave(plot = integ_plot, filename = file.path(output_dir,"Integrated_split_by_origidentx.pdf"), width = 6, height = 3)

# ==============================================================================
# Figure 3B: EC-subcluster marker expression
# ==============================================================================

celltype.names <- levels(integ_data)
tcell.names <- grep("EC", celltype.names,value = TRUE)
tcells_integ <- subset(integ_data, idents = tcell.names)

EC_cols <- c('EC cells-5'='#ABA300',
             'EC cells-4'='#E68613',
             'EC cells-1'='#F8766D',
             'EC cells-3'='#00BE67',
             'EC cells-2'='#00BCD7'
             )

EC_cell_marker<- VlnPlot(tcells_integ, features = c("LMX1A", "MNX1", "MAFB", "EBF1"), ncol = 2, pt.size=0,assay = 'SCT', cols = EC_cols)
ggsave(plot = EC_cell_marker,filename =  file.path(output_dir,"EC_marker_violinplot_ECclusters.jpg"))

# ==============================================================================
# Figure 3C: TF expression and motif activity
# ==============================================================================

DefaultAssay(integ_data) <- 'ATAC'
motif.name <- ConvertMotifID(integ_data, name = 'EBF1')
gene_plot_ebf1 <- FeaturePlot(integ_data, features = "sct_EBF1", reduction = 'umap',
                              min.cutoff = 'q10',max.cutoff = 'q90') + NoLegend() + NoAxes() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
motif_plot_ebf1 <- FeaturePlot(integ_data, features = motif.name, cols = c("lightgrey", "darkred"), reduction = 'umap',
                        min.cutoff = 'q10',max.cutoff = 'q90')+ NoLegend() + NoAxes() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
motif.name <- ConvertMotifID(integ_data, name = 'LMX1A')
gene_plot_lmx1a <- FeaturePlot(integ_data, features = "sct_LMX1A", reduction = 'umap',
                              min.cutoff = 'q10',max.cutoff = 'q90') + NoLegend() + NoAxes() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
motif_plot_lmx1a <- FeaturePlot(integ_data, features = motif.name, cols = c("lightgrey", "darkred"), reduction = 'umap',
                        min.cutoff = 'q10',max.cutoff = 'q90')+ NoLegend() + NoAxes() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
motif.name <- ConvertMotifID(integ_data, name = 'MNX1')
gene_plot_MNX1<- FeaturePlot(integ_data, features = "sct_MNX1", reduction = 'umap',
                             min.cutoff = 'q10',max.cutoff = 'q90')+ NoLegend() + NoAxes() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
motif_plot_MNX1 <- FeaturePlot(integ_data, features = motif.name, cols = c("lightgrey", "darkred"), reduction = 'umap',
                        min.cutoff = 'q10',max.cutoff = 'q90') + NoLegend() + NoAxes() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

ggsave(filename = file.path(output_dir,"LMX1A_EBF1_MNX1_TFs_TopTFs_EC_branch.pdf"),width = 6,height = 12,
                plot = gridExtra::marrangeGrob(list(gene_plot_ebf1, gene_plot_lmx1a, gene_plot_MNX1, motif_plot_ebf1, motif_plot_lmx1a, motif_plot_MNX1), nrow = 3, ncol = 2), device = "pdf")
motif.name <- ConvertMotifID(integ_data, name = 'EBF1')
motif_EBF1<- MotifPlot(
  object = integ_data,
  motifs = motif.name)
ggsave(motif_EBF1, filename =file.path(output_dir,"motif_EBF1_logo.png"),width = 2.14, height = 1.51, dpi = 300)

motif.name <- ConvertMotifID(integ_data, name = 'LMX1A')
motif_LMX1A<- MotifPlot(
  object = integ_data,
  motifs = motif.name)
ggsave(motif_LMX1A, filename =file.path(output_dir,"motif_LMX1A_logo.png"),width = 2.14, height = 1.51, dpi = 300)

motif.name <- ConvertMotifID(integ_data, name = 'MNX1')
motif_MNX1 <- MotifPlot(
  object = integ_data,
  motifs = motif.name)
ggsave(motif_MNX1, filename =file.path(output_dir,"motif_MNX1_logo.png"),width = 2.14, height = 1.51, dpi = 300)
           
# ==============================================================================
# Figure 3D: Average expression of cell-type marker genes
# ==============================================================================

Idents(integ_data) <- 'celltype2'
cell.list <- WhichCells(integ_data, idents = c("Beta cells","Hi INS Beta cells", "Alpha cells", "SST cells", "EC cells-1", "EC cells-2","EC cells-3", "EC cells-4","EC cells-5"))
integ_subset <- integ_data[, cell.list]
sig_markers <- c( "EBF1","ERO1B","SLC30A8", "ROBO1","ROBO2","HCN1","CADM1", "G6PC2", "PLAGL1","INS",  "IAPP","MAFB", "FTL", "TTR",
                "GCG", "NRXN3","UNC5C","DPP4", "SST", "PCDH9",
                "LMX1A", "LMX1B","GRIA1","PPFIA2", "HDAC9", "AFF3", "TPH1", "ZBTB7C", "GRID2","LRP1B", "SLIT1")
integ_subset<- NormalizeData(object = integ_subset, assay = "RNA")                
avg <- AverageExpression(integ_subset,group.by = "celltype2",features = sig_markers,assays = "RNA")               
avg_scaled = t(scale(t(avg$RNA))) 
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
marker_heatmap <- ComplexHeatmap::Heatmap(
  avg_scaled,
  heatmap_legend_param = list(
    title = "Normalized gene expression",
    title_position = "leftcenter-rot"
  ),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = TRUE,
  show_row_names = TRUE,
  col = col_fun,
  column_title = "Cell-type-specific marker average expression"
)

pdf(
  file.path(output_dir, "average_marker_expression.pdf"),
  width = 8,
  height = 8
)

ComplexHeatmap::draw(marker_heatmap)

dev.off()
              
# ==============================================================================
# Figure 3E: High-INS beta-cell marker expression
# ==============================================================================

hiINS_marker_expressions<- VlnPlot(integ_data, features = c("INS", "eGFP", "IAPP", "ISL1"), ncol = 2, pt.size=0,assay = 'SCT', cols = my_cols_int)
ggsave(plot = hiINS_marker_expressions,filename =  file.path(output_dir,"HiINS_marker_violinpplot_onintegrated.jpg"), width = 10,height = 8)         






       