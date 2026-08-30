#!/bin/bash

#SBATCH --job-name=pyscenic
#SBATCH --account=st-wasser-1
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=256G
#SBATCH --output=pyscenic_%j.out
#SBATCH --error=pyscenic_%j.err

# ==============================================================================
# SCENIC analysis: Regulon inference and activity
#
# Description:
# Performs pySCENIC analysis of integrated EC and beta-cell populations.
# The Seurat object is converted to loom format, GRNBoost2 and cisTarget are
# run 20 times, and TF-target interactions identified in at least 80% of runs
# are retained as consensus regulons. AUCell is then used to quantify
# regulon activity across cells.
#
# Input:
# data/processed/ECandBetasubs_integrated_seurat.rds
#
# Reference files:
# hs_hgnc_tfs.txt
# hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
# motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
#
# Outputs:
# ec_and_beta_integ.loom
# regulons.csv
# ec_and_beta_integ_AUCell.loom
#
# Environment:
# R package versions are recorded in renv.lock.
# pySCENIC was run using the pyscenic_test conda environment.
# ==============================================================================


# ==============================================================================
# 1. Prepare loom input
# ==============================================================================

Rscript - <<'EOF'

library(Seurat)
library(SeuratDisk)

ECandBetasubs <- readRDS(
    "data/processed/ECandBetasubs_integrated_seurat.rds"
)

DefaultAssay(ECandBetasubs) <- "RNA"

loom <- SeuratDisk::as.loom(
    ECandBetasubs,
    filename = "ec_and_beta_integ.loom"
)

loom$close_all()

EOF


# ==============================================================================
# 2. GRN inference and motif enrichment
# ==============================================================================

source ~/.bashrc
conda activate pyscenic_test

subset="ec_and_beta_integ"

list_TFs="hs_hgnc_tfs.txt"
ranking_db="hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
motif_ann="motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"

for repeat in {1..20}
do
    echo "Processing pySCENIC repeat ${repeat}"

    pyscenic grn \
        "${subset}.loom" \
        "${list_TFs}" \
        --output "${subset}_adjacencies_${repeat}.csv" \
        --num_workers 4 \
        --seed "${repeat}"

    pyscenic ctx \
        "${subset}_adjacencies_${repeat}.csv" \
        "${ranking_db}" \
        --annotations_fname "${motif_ann}" \
        --expression_mtx_fname "${subset}.loom" \
        --output "${subset}_regulons_${repeat}.csv" \
        --num_workers 4
done


# ==============================================================================
# 3. Consensus regulon selection
# ==============================================================================

# TF-target interactions identified in at least 80% (>=16/20) of the
# independent pySCENIC runs were retained and saved as:
#
# regulons.csv


# ==============================================================================
# 4. Regulon activity scoring
# ==============================================================================

pyscenic aucell \
    "${subset}.loom" \
    "regulons.csv" \
    --output "${subset}_AUCell.loom" \
    --num_workers 4

conda deactivate