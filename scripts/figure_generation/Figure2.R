# ==============================================================================
# Figure 2: Characterization of cell types at Day 30 and Day 50
#
# Description:
# Generates cell-type proportion plots, cell-type marker heatmaps,
# beta-cell GO enrichment plots, and differential gene expression and
# motif activity plots for Day 30 and Day 50 samples.
#
# Inputs:
#   Processed and cell-type-annotated Day 30 multiome dataset
#   Processed and cell-type-annotated Day 50 multiome dataset

# Outputs:
#   figures/Figure2/
#
# Environment:
#   R package versions are recorded in renv.lock.
#
# NOTE:
#   Figure 2F/G volcano analyses require additional review of significance
#   filtering and feature-label selection before final manuscript submission.
# ============================================================================
# Packages and settings
# ==============================================================================
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(EnhancedVolcano)
set.seed(1234)

# ==============================================================================
# Input data and output directory
# ==============================================================================

d30_input <- readRDS(
  "data/processed/day30_annotated_multiome.rds"
)

d50_input <- readRDS(
  "data/processed/day50_annotated_multiome.rds"
)


output_dir <- "figures/Figure2"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Shared plotting settings
# ==============================================================================

celltype_colors <- c(
  "Beta cells" = "#FF61CC",
  "Alpha cells" = "#00C19A",
  "EC cells" = "#E68613",
  "Polyhormonal cells" = "#8494FF",
  "SST cells" = "#0CB702",
  "Alpha-Beta cells" = "#00A9FF",
  "UNK" = "#D3D3D3"
)

# ==============================================================================
# Figure 2A: Cell-type proportions at Day 30 and Day 50
# ==============================================================================
# Consolidate Day 30 cell types
d30_prop_obj <- RenameIdents(
  d30_input,
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

# Consolidate Day 50 cell types
d50_prop_obj <- RenameIdents(
  d50_input,
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

# Calculate percentages
df30 <- as.data.frame(prop.table(table(Idents(d30_prop_obj))) * 100)
colnames(df30) <- c("Clusters", "Proportions")
df30$Dataset <- "Day30"

df50 <- as.data.frame(prop.table(table(Idents(d50_prop_obj))) * 100)
colnames(df50) <- c("Clusters", "Proportions")
df50$Dataset <- "Day50"

# Combine datasets
data <- rbind(df30, df50)

data$Clusters <- factor(
  data$Clusters,
  levels = c(
    "Beta cells",
    "Alpha cells",
    "EC cells",
    "Polyhormonal cells",
    "SST cells",
    "Alpha-Beta cells",
    "UNK"
  )
bp<- ggplot(data, aes(x=Dataset, y=Proportions, fill=Clusters))+ geom_bar(width = 0.8, stat = "identity")+ theme_minimal() + scale_fill_manual(values=celltype_colors)

pdf(file.path(output_dir, "celltype_proportions.pdf"), width = 4, height = 4)
print(bp)
dev.off()


# ==============================================================================
# Figure 2B: Day 30 cell-type marker heatmap
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
d30_heatmap_markers_4majorclusters <- FindAllMarkers(d30_input_data2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay="RNA")
#write.csv(d30_heatmap_markers_4majorclusters, "d30_heatmap_markers_4majorclusters.csv")
d30_markers_backup <- d30_heatmap_markers_4majorclusters

d30_markers_backup %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>% arrange(desc(avg_log2FC), .by_group = TRUE)%>%
    slice_head(n = 10) %>%
    ungroup() -> top10
top10_sub<- top10[!top10$cluster %in% c("Alpha-Beta cells", "UNK", "Polyhormonal cells"), ]

d30_input_data2_subset_forplot <- subset(d30_input_data2, idents = c("Beta cells", "Alpha cells", "EC cells", "SST cells"))
day30_heatmap<- DoHeatmap(d30_input_data2_subset_forplot, features = top10_sub$gene,group.colors = c("#FF61CC", '#00C19A','#E68613','#0CB702')) + scale_fill_gradientn(colors = c("blue", "white", "red")) 

pdf(file.path(output_dir, "heatmap_d30_markers_top10_10x20.pdf"),width = 5, height = 8)
print(day30_heatmap)
dev.off()

# ==============================================================================
# Figure 2C: Day 50 cell-type marker heatmap
# ==============================================================================
d50_input_heatmap <- RenameIdents(object=d50_input, 
                           "Beta cells-1" = "Beta cells",
                           "Hi INS Beta cells-2" = "Beta cells",
                           "Alpha-Beta cells"="Alpha-Beta cells",
                           "Alpha cells"= "Alpha cells",
                           "SST+ cells" = "SST cells", 
                           "Polyhormonal cells"="Polyhormonal cells",
                           "EC cells-1"="EC cells",
                           "EC cells-2"="EC cells",
                           "EC cells-3"="EC cells",
                           "EC cells-4"="EC cells",
                           "UNK"= "UNK")
tcells_integ_heatmap_D50 <- subset(d50_input_heatmap, idents = c("Beta cells", "Alpha cells", "EC cells", "SST cells"))

D50_heatmap_markers<- FindAllMarkers(tcells_integ_heatmap_D50, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay="RNA")
#write.csv(D50_heatmap_markers,"d50_markers_EC_beta_subclusters_merged.csv")
DefaultAssay(tcells_integ_heatmap_D50)<- "RNA"
D50_heatmap_markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>% arrange(desc(avg_log2FC), .by_group = TRUE)%>%
    slice_head(n = 10) %>%
    ungroup() -> top10_50

day50_heatmap<- DoHeatmap(tcells_integ_heatmap_D50, features = top10_50$gene,group.colors = c("#FF61CC", '#00C19A','#E68613','#0CB702')) + scale_fill_gradientn(colors = c("blue", "white", "red")) 

pdf(file.path(output_dir, "Day50_markers_heatmap_top10.pdf"), width = 8, height = 10)
print(day50_heatmap)
dev.off()

# ==============================================================================
# Figure 2D: Day 30 beta-cell pathway enrichment
# ==============================================================================
d30_inputplot <- RenameIdents(object=d30_input, 
                           "Beta cells-1" = "Beta cells",
                           "Hi INS Beta cells-2" = "Beta cells",
                           "Alpha-Beta cells"="Alpha-Beta cells",
                           "Alpha cells"= "Alpha cells",
                           "SST+ cells" = "SST cells", 
                           "Polyhormonal cells"="Polyhormonal cells",
                           "EC cells-1"="EC cells-1",
                           "EC cells-2"="EC cells-2",
                           "EC cells-3"="EC cells-3",
                           "EC cells-4"="EC cells-4",
                           "UNK"= "UNK")
DefaultAssay(d30_inputplot) <- "RNA"
d30_beta_markers <- FindMarkers(d30_inputplot, ident.1 = "Beta cells")                           
d30_beta_markers<- d30_beta_markers %>% arrange(desc(avg_log2FC))
D30_beta_topgenesup <- d30_beta_markers[d30_beta_markers$avg_log2FC > 0,]
d30_beta_topgenenames<- rownames(D30_beta_topgenesup)                           
topgenes_beta_d30 <- enrichGO(gene = d30_beta_topgenenames, ont = "BP",
                   OrgDb =org.Hs.eg.db,keyType= 'SYMBOL',
                   pvalueCutoff = 0.05)                           
    
d30_selected_pathways <- c(
  "axonogenesis",
  "axon development",
  "modulation of chemical synaptic transmission",
  "regulation of nervous system development",
  "regulation of protein secretion",
  "response to endoplasmic reticulum stress",
  "peptide hormone secretion",
  "peptide transport"
)

d30_upregulated_paths_beta_plot <- topgenes_beta_d30

d30_upregulated_paths_beta_plot@result <-
  topgenes_beta_d30@result[
    topgenes_beta_d30@result$Description %in% d30_selected_pathways,
  ]

path_d30_up<- barplot(d30_upregulated_paths_beta_plot,)  +ggtitle('Enriched Gene Sets in Day 30 Beta cells ')
ggsave(path_d30_up,filename = file.path(
    output_dir,"D30_Beta_cluster_upregulated_pathways.jpg"),height = 4, width = 6)

# ==============================================================================
# Figure 2E: Day 50 beta-cell pathway enrichment
# ==============================================================================
d50_inputplot <- RenameIdents(object=d50_input, 
                           "Beta cells-1" = "Beta cells",
                           "Hi INS Beta cells-2" = "Beta cells",
                           "Alpha-Beta cells"="Alpha-Beta cells",
                           "Alpha cells"= "Alpha cells",
                           "SST+ cells" = "SST cells", 
                           "Polyhormonal cells"="Polyhormonal cells",
                           "EC cells-1"="EC cells-1",
                           "EC cells-2"="EC cells-2",
                           "EC cells-3"="EC cells-3",
                           "EC cells-4"="EC cells-4",
                           "UNK"= "UNK")

DefaultAssay(d50_inputplot) <- "RNA"

d50_beta_markers <- FindMarkers(
  d50_inputplot,
  ident.1 = "Beta cells"
)

d50_beta_markers <- d50_beta_markers %>%
  arrange(desc(avg_log2FC))
# NOTE: Day 50 beta markers were filtered at adjusted p < 0.005
# before GO enrichment, consistent with the original analysis.
d50_beta_markers_0005 <- d50_beta_markers[
  d50_beta_markers$p_val_adj < 0.005,
]

D50_beta_topgenesup <- d50_beta_markers_0005[
  d50_beta_markers_0005$avg_log2FC > 0,
]

D50_beta_topgenenames <- rownames(D50_beta_topgenesup)

topgenes_beta_Day50 <- enrichGO(
  gene = D50_beta_topgenenames,
  ont = "BP",
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  pvalueCutoff = 0.05
)             

d50_selected_pathways <- c(
  "axonogenesis",
  "protein folding in endoplasmic reticulum",
  "response to endoplasmic reticulum stress",
  "regulation of protein secretion",
  "modulation of chemical synaptic transmission",
  "protein localization to extracellular region",
  "regulation of nervous system development",
  "regulation of insulin secretion"
)

d50_upregulated_paths_beta_plot <- topgenes_beta_Day50

d50_upregulated_paths_beta_plot@result <-
  topgenes_beta_Day50@result[
    topgenes_beta_Day50@result$Description %in% d50_selected_pathways,
  ]

path_d50_up <- barplot(
  d50_upregulated_paths_beta_plot
) +
  ggtitle("Enriched Gene Sets in Day 50 Beta cells")

ggsave(path_d50_up,filename = file.path(
    output_dir,"D50_Beta_cluster_upregulated_pathways.jpg"),height = 4, width = 6)
                  
                   
# ==============================================================================
# Figure 2F: Day 30 differential expression and motif activity
# ==============================================================================

# NOTE: Statistical filtering and volcano labeling to be reviewed.

DefaultAssay(d30_input)<- "RNA"
D30_alpha_vs_beta <- FindMarkers(d30_input, ident.1 = c("Beta cells-1", "Hi INS Beta cells-2"), 
                                 ident.2 = "Alpha cells",
                                 only.pos = FALSE, test.use = 'LR')
D30_alpha_vs_beta_sig<- D30_alpha_vs_beta[D30_alpha_vs_beta$p_val<0.05, ]
D30_alpha_vs_beta_sig<- D30_alpha_vs_beta_sig %>% arrange(desc(avg_log2FC))
beta_names<- rownames(D30_alpha_vs_beta_sig[D30_alpha_vs_beta_sig$avg_log2FC>1.5, ])
alpha_names<- rownames(D30_alpha_vs_beta_sig[D30_alpha_vs_beta_sig$avg_log2FC< -1.2, ])
labs<- c(beta_names, alpha_names)
add<- c("INS", "eGFP")
selectlab<- c(labs,add)

d30volcano_alpha_beta <- EnhancedVolcano(D30_alpha_vs_beta_sig,
    lab = rownames(D30_alpha_vs_beta_sig),drawConnectors = TRUE,
    widthConnectors = 0.75,
    x = 'avg_log2FC',selectLab= selectlab,
    y = 'p_val_adj',pointSize = 3.0,labSize = 4.0, title='Alpha cells vs Beta cells D30')
d30volcano_alpha_beta

pdf("D30_Volcanoplot_DEGs_Alpha_vs_Beta2.pdf", width=10, height=10)
print(d30volcano_alpha_beta)
dev.off()

DefaultAssay(d30_input)<- 'chromvar'
differential.activity <- FindMarkers(
  object = d30_input,
    ident.1 = c("Beta cells-1", "Hi INS Beta cells-2"),
    ident.2 = "Alpha cells", 
    only.pos = FALSE,
    mean.fxn = rowMeans,
    fc.name = "avg_diff"
)
## Day 30 -DAR volcano
DefaultAssay(d30_input)<- "ATAC"
differential.activity$gene <- ConvertMotifID(d30_input, id = rownames(differential.activity))

differential.activity_sig <- differential.activity[differential.activity$p_val<0.05, ]
differential.activity_sig<- differential.activity_sig %>% arrange(desc(avg_diff))
d30volcano_alpha_beta_motif <- EnhancedVolcano(differential.activity_sig,
    lab = differential.activity_sig$gene, selectLab= selectlab,
    x = 'avg_diff',FCcutoff = 0.5,drawConnectors = TRUE,widthConnectors = 0.75,
    y = 'p_val_adj',pointSize = 2.0,labSize = 5, title='Differential Motif Activity Beta cells vs Alpha cells Day 30')
d30volcano_alpha_beta_motif

pdf("D30_Volcanoplot_Differential_motif_activity_Alpha_vs_Beta.pdf", width=10, height=10)
print(d30volcano_alpha_beta_motif)
dev.off()                   
                   
# ==============================================================================
# Figure 2G: Day 50 differential expression and motif activity
# ==============================================================================

# NOTE: Statistical filtering and volcano labeling to be reviewed.                 
                   
DefaultAssay(d50_input)<- "RNA"
D50_alpha_vs_beta <- FindMarkers(d50_input, ident.1 = c("Beta cells-1", "Hi INS Beta cells-2"), 
                                 ident.2 = "Alpha cells",
                                 only.pos = FALSE, test.use = 'LR')                   
                   
d50_alpha_vs_beta_sig<- D50_alpha_vs_beta[D50_alpha_vs_beta$p_val<0.05, ]
d50_alpha_vs_beta_sig<- d50_alpha_vs_beta_sig %>% arrange(desc(avg_log2FC))                  
beta_names50<- rownames(D30_alpha_vs_beta_sig[D30_alpha_vs_beta_sig$avg_log2FC>1.8, ])
alpha_names50<- rownames(D30_alpha_vs_beta_sig[D30_alpha_vs_beta_sig$avg_log2FC< -1.4, ])

length(beta_names50)
length(alpha_names50)

labs<- c(beta_names50, alpha_names50)
add<- c("INS", "eGFP")
selectlab<- c(labs,add)
                   
d50volcano_alpha_beta <- EnhancedVolcano(d50_alpha_vs_beta_sig,
    lab = rownames(d50_alpha_vs_beta_sig),
    x = 'avg_log2FC', selectLab= selectlab,drawConnectors = TRUE,
    widthConnectors = 0.75,
    y = 'p_val_adj',pointSize = 3.0,labSize = 4.0,title='Alpha cells vs Beta cells D50' )
d50volcano_alpha_beta


pdf("D50_Volcanoplot_DEGs_Alpha_vs_Beta2.pdf", width=10, height=10)
print(d50volcano_alpha_beta)
dev.off() 

## Day 50 DAR 
DefaultAssay(d50_input)<- 'chromvar'
differential.activity_d50 <- FindMarkers(
  object = d50_input,
    ident.1 = c("Beta cells-1", "Hi INS Beta cells-2"),
    ident.2 = "Alpha cells", 
    only.pos = FALSE,
    mean.fxn = rowMeans,
    fc.name = "avg_diff"
)

DefaultAssay(d50_input)<- "ATAC"
differential.activity_d50$gene <- ConvertMotifID(d50_input, id = rownames(differential.activity_d50))

differential.activityD50_sig <- differential.activity_d50[differential.activity_d50$p_val<0.05, ]
differential.activityD50_sig<- differential.activityD50_sig %>% arrange(desc(avg_diff))     
D50volcano_alpha_beta_motif <- EnhancedVolcano(differential.activityD50_sig,
    lab = differential.activityD50_sig$gene, selectLab= selectlab,
    x = 'avg_diff',FCcutoff = 0.5,drawConnectors = TRUE,widthConnectors = 0.75,
    y = 'p_val_adj',pointSize = 2.0,labSize = 5, title='Differential Motif Activity Beta cells vs Alpha cells Day 50')
D50volcano_alpha_beta_motif

pdf("D50_Volcanoplot_Differential_motif_activity_Alpha_vs_Beta.pdf", width=10, height=10)
print(D50volcano_alpha_beta_motif)
dev.off()    


                    