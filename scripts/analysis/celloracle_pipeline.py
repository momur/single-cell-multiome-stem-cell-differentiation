# ==============================================================================
# CellOracle analysis: GRN inference and in silico TF perturbation
#
# Description:
# Constructs cell-type-specific gene regulatory networks for integrated EC and
# beta-cell populations using CellOracle and performs in silico transcription
# factor knockout simulations.
#
# Inputs:
# data/processed/EC_BETA_integ_celloracle.h5ad
# data/processed/base_GRN_dataframe.parquet
#
# Outputs:
# Results used for Figure 4E-F and associated supplementary analyses.

#
# Environment:
# Python package versions are recorded in the project environment.
# ==============================================================================


# ==============================================================================
# Packages and settings
# ==============================================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import celloracle as co

np.random.seed(1234)


# ==============================================================================
# Input data and output directories
# ==============================================================================

adata_path = "data/processed/EC_BETA_integ_celloracle.h5ad"
base_grn_path = "data/processed/base_GRN_dataframe.parquet"

processed_dir = "data/processed"
figure_dir = "figures/Figure4"

os.makedirs(processed_dir, exist_ok=True)
os.makedirs(figure_dir, exist_ok=True)


# ==============================================================================
# Load integrated EC and beta-cell dataset
# ==============================================================================

adata = sc.read_h5ad(adata_path)

print(f"Number of cells: {adata.shape[0]}")
print(f"Number of genes: {adata.shape[1]}")

# Preserve raw counts for CellOracle analysis
adata.layers["counts"] = adata.X.copy()


# ==============================================================================
# Initialize CellOracle
# ==============================================================================

oracle = co.Oracle()

oracle.import_anndata_as_raw_count(
    adata=adata,
    cluster_column_name="celltype",
    embedding_name="X_umap"
)


# ==============================================================================
# Import ATAC-informed base GRN
# ==============================================================================

# Load custom ATAC-informed base GRN generated from our multiome dataset.
# The base GRN was constructed using Cicero co-accessibility, hg38 TSS
# annotation, and TF motif scanning.

base_GRN = pd.read_parquet(
    base_grn_path,
    engine="pyarrow"
)

oracle.import_TF_data(
    TF_info_matrix=base_GRN
)


# ==============================================================================
# PCA and KNN imputation
# ==============================================================================

oracle.perform_PCA()

# Select the number of PCs using the explained-variance profile,
# with a maximum of 50 PCs.
n_comps = np.where(
    np.diff(
        np.diff(
            np.cumsum(oracle.pca.explained_variance_ratio_)
        ) > 0.002
    )
)[0][0]

n_comps = min(n_comps, 50)

# Number of neighbors corresponding to 2.5% of cells
n_cells = oracle.adata.shape[0]
k = int(0.025 * n_cells)

oracle.knn_imputation(
    n_pca_dims=n_comps,
    k=k,
    balanced=True,
    b_sight=k * 8,
    b_maxl=k * 4,
    n_jobs=4
)

oracle.to_hdf5(
    os.path.join(
        processed_dir,
        "EC_beta_subset.celloracle.oracle"
    )
)


# ==============================================================================
# Cell-type-specific GRN inference
# ==============================================================================

links = oracle.get_links(
    cluster_name_for_GRN_unit="celltype",
    alpha=10,
    verbose_level=10
)

links.filter_links(
    p=0.001,
    weight="coef_abs",
    threshold_number=2000
)

links.get_network_score()

links.to_hdf5(
    file_path=os.path.join(
        processed_dir,
        "links.celloracle.links"
    )
)


# ==============================================================================
# Export filtered cell-type-specific GRNs
# ==============================================================================

cell_types = [
    "EC cells-5",
    "EC cells-4",
    "EC cells-3",
    "EC cells-2",
    "EC cells-1",
    "Beta cells-1",
    "Hi INS Beta cells-2"
]

grn_output_dir = os.path.join(
    processed_dir,
    "celltype_GRNs"
)

os.makedirs(grn_output_dir, exist_ok=True)

for cell_type in cell_types:

    output_name = cell_type.replace(" ", "_")

    links.filtered_links[cell_type].to_csv(
        os.path.join(
            grn_output_dir,
            f"filtered_GRN_{output_name}.csv"
        ),
        index=False
    )


# ==============================================================================
# Prepare GRNs for in silico perturbation
# ==============================================================================

links.filter_links(
    p=0.001,
    weight="coef_abs",
    threshold_number=2000
)

oracle.get_cluster_specific_TFdict_from_Links(
    links_object=links
)

oracle.fit_GRN_for_simulation(
    alpha=10,
    use_cluster_specific_TFdict=True
)


# ==============================================================================
# In silico TF knockout simulations
# ==============================================================================

# TFs evaluated in the perturbation analysis
tf_list = [
    "ZBTB7C",
    "FOXA2",
    "EBF1",
    "ESR1",
    "STAT1",
    "PLAGL1",
    "LMX1A",
    "FOXO1"
]

simulation_dir = os.path.join(
    figure_dir,
    "CellOracle_perturbations"
)

os.makedirs(simulation_dir, exist_ok=True)


for goi in tf_list:

    print(f"Running CellOracle simulation for {goi} KO")

    # ----------------------------------------------------------
    # Simulate TF knockout
    # ----------------------------------------------------------

    oracle.simulate_shift(
        perturb_condition={goi: 0.0},
        n_propagation=3
    )

    # ----------------------------------------------------------
    # Estimate cell-state transition probabilities
    # ----------------------------------------------------------

    oracle.estimate_transition_prob(
        n_neighbors=200,
        knn_random=True,
        sampled_fraction=1
    )

    oracle.calculate_embedding_shift(
        sigma_corr=0.05
    )

    # ----------------------------------------------------------
    # Calculate simulation flow on grid
    # ----------------------------------------------------------

    n_grid = 40

    oracle.calculate_p_mass(
        smooth=0.8,
        n_grid=n_grid,
        n_neighbors=200
    )

    oracle.calculate_mass_filter(
        min_mass=10,
        plot=False
    )

    # ----------------------------------------------------------
    # Plot simulated cell-state shifts
    # ----------------------------------------------------------

    fig, ax = plt.subplots(figsize=(8, 8))

    oracle.plot_cluster_whole(
        ax=ax,
        s=10
    )

    oracle.plot_simulation_flow_on_grid(
        scale=30,
        ax=ax,
        show_background=False
    )

    ax.set_title(
        f"Simulated cell identity shift: {goi} KO"
    )

    plt.savefig(
        os.path.join(
            simulation_dir,
            f"{goi}_KO_simulation.png"
        ),
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()


# ==============================================================================
# Save final CellOracle object
# ==============================================================================

oracle.to_hdf5(
    os.path.join(
        processed_dir,
        "EC_beta_subset_simulations.celloracle.oracle"
    )
)