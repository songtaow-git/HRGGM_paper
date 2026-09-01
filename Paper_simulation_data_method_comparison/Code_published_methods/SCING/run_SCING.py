"""
SCING consensus-network inference.

The input matrix is interpreted as local features by samples. SCING builds
multiple subsampled directed networks, aggregates edge weights across runs,
and saves directed and symmetrized consensus weight matrices.
"""

from __future__ import annotations

import os
import shutil
import sys
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd


# =============================================================================
# Analysis settings
# =============================================================================

# All relative paths are resolved from the directory containing this script.
SCRIPT_DIR = Path(__file__).resolve().parent

# Directory containing the SCING source tree.
SCING_REPO_DIR = SCRIPT_DIR / "SCING-main"

# Numeric input matrix with features in rows and samples in columns.
INPUT_MATRIX_FILE = SCRIPT_DIR / "Data_whole_rep1.txt"

# Final output directory.
OUTPUT_DIR = SCRIPT_DIR / "SCING_Output_rep1"

# Number of independently subsampled SCING networks.
N_NETWORKS = 100

# Optional pre-subsampling limit. Use None to retain all samples.
MAX_CELLS_FOR_SCING: int | None = None

# Minimum fraction of intermediate networks in which an edge must appear.
MINIMUM_EDGE_APPEARANCE_THRESHOLD = 0.0

# Fraction of samples used by each intermediate SCING network.
SUBSAMPLE_PERC = 0.7

# SCING neighborhood size and number of principal components.
NNEIGHBORS = 100
NPCS = 20

# Computational resources passed to SCING.
NCORE = 1
MEM_PER_CORE = int(16e9)

# Random seed used for optional pre-subsampling and reproducible execution.
RANDOM_SEED = 1

# Method used to combine opposite directed weights into one undirected weight.
# Supported values are "max", "mean", and "sum".
UNDIRECTED_COMBINE_METHOD = "max"

# Final output files.
DIRECTED_WEIGHT_MATRIX_FILE = OUTPUT_DIR / "SCING_weight_matrix.txt"
UNDIRECTED_WEIGHT_MATRIX_FILE = (
    OUTPUT_DIR / "SCING_weight_matrix_undirected_max.txt"
)


# =============================================================================
# SCING import
# =============================================================================

def add_scing_to_path() -> None:
    """Add the SCING source directory to the Python import path."""
    src_dir = SCING_REPO_DIR / "src"
    if not src_dir.is_dir():
        raise FileNotFoundError(f"SCING source directory not found: {src_dir}")

    src_text = str(src_dir)
    if src_text not in sys.path:
        sys.path.insert(0, src_text)


add_scing_to_path()

try:
    import anndata as ad
    from scing import build
except Exception as exc:
    raise ImportError(
        "SCING and its Python dependencies must be installed before running "
        "this script."
    ) from exc


# =============================================================================
# Input and matrix utilities
# =============================================================================

def read_numeric_matrix(filename: Path) -> np.ndarray:
    """Read a numeric matrix from a whitespace- or comma-delimited text file."""
    if not filename.is_file():
        raise FileNotFoundError(f"Input matrix not found: {filename}")

    read_attempts = (
        {"sep": r"\s+", "engine": "python"},
        {"sep": ",", "engine": "c"},
    )

    for options in read_attempts:
        try:
            matrix = pd.read_csv(
                filename,
                header=None,
                **options,
            ).to_numpy(dtype=float)
            if matrix.size > 0:
                break
        except (ValueError, pd.errors.ParserError):
            matrix = np.empty((0, 0), dtype=float)
    else:
        matrix = np.empty((0, 0), dtype=float)

    if matrix.size == 0:
        raise ValueError(f"No numeric data were read from: {filename}")

    if not np.isfinite(matrix).all():
        raise ValueError("The input matrix contains NaN or infinite values.")

    return matrix


def subsample_samples(
    matrix: np.ndarray,
    maximum_samples: int | None,
    random_seed: int,
) -> np.ndarray:
    """Optionally retain a random subset of sample columns."""
    if maximum_samples is None or matrix.shape[1] <= maximum_samples:
        return matrix

    if maximum_samples < 1:
        raise ValueError("MAX_CELLS_FOR_SCING must be positive or None.")

    rng = np.random.default_rng(random_seed)
    selected = np.sort(
        rng.choice(matrix.shape[1], size=maximum_samples, replace=False)
    )
    return matrix[:, selected]


def validate_feature_variance(matrix: np.ndarray) -> None:
    """Require finite, nonzero variance for every feature."""
    feature_sd = np.std(matrix, axis=1)
    invalid = np.flatnonzero(~np.isfinite(feature_sd) | (feature_sd <= 1e-12))

    if invalid.size > 0:
        raise ValueError(
            "The input contains constant or invalid feature rows. "
            f"First affected row: {invalid[0] + 1}"
        )


def generate_feature_names(number_of_features: int) -> np.ndarray:
    """Generate feature labels used by SCING intermediate edge tables."""
    return np.asarray(
        [f"Gene_{index + 1}" for index in range(number_of_features)],
        dtype=str,
    )


# =============================================================================
# Intermediate SCING networks
# =============================================================================

def standardize_edge_columns(edge_table: pd.DataFrame) -> pd.DataFrame:
    """Map a SCING edge table to source, target, and importance columns."""
    table = edge_table.copy()
    unnamed = [
        column
        for column in table.columns
        if str(column).startswith("Unnamed")
    ]
    if unnamed:
        table = table.drop(columns=unnamed)

    lower_map = {str(column).lower(): column for column in table.columns}

    source_target_candidates = (
        ("source", "target"),
        ("regulator", "target"),
        ("tf", "target"),
        ("from", "to"),
    )

    source_column = None
    target_column = None

    for source_name, target_name in source_target_candidates:
        if source_name in lower_map and target_name in lower_map:
            source_column = lower_map[source_name]
            target_column = lower_map[target_name]
            break

    if source_column is None or target_column is None:
        gene_like_counts = []
        for column in table.columns:
            values = table[column].astype(str)
            count = values.str.match(r"^Gene_\d+$").sum()
            gene_like_counts.append((column, count))

        gene_like_counts.sort(key=lambda item: item[1], reverse=True)

        if len(gene_like_counts) < 2 or gene_like_counts[1][1] == 0:
            raise ValueError(
                "Source and target columns could not be identified in a "
                "SCING intermediate edge table."
            )

        source_column = gene_like_counts[0][0]
        target_column = gene_like_counts[1][0]

    importance_column = None
    for candidate in ("importance", "weight", "score", "importance_x"):
        if candidate in lower_map:
            importance_column = lower_map[candidate]
            break

    if importance_column is None:
        numeric_candidates = []
        for column in table.columns:
            if column in (source_column, target_column):
                continue
            numeric_count = pd.to_numeric(
                table[column],
                errors="coerce",
            ).notna().sum()
            if numeric_count > 0:
                numeric_candidates.append((column, numeric_count))

        if not numeric_candidates:
            raise ValueError(
                "An edge-weight column could not be identified in a SCING "
                "intermediate edge table."
            )

        numeric_candidates.sort(key=lambda item: item[1], reverse=True)
        importance_column = numeric_candidates[0][0]

    standardized = table[
        [source_column, target_column, importance_column]
    ].copy()
    standardized.columns = ["source", "target", "importance"]

    standardized["source"] = standardized["source"].astype(str)
    standardized["target"] = standardized["target"].astype(str)
    standardized["importance"] = pd.to_numeric(
        standardized["importance"],
        errors="coerce",
    )

    standardized = standardized.dropna(
        subset=["source", "target", "importance"]
    )
    standardized = standardized[
        standardized["source"].str.match(r"^Gene_\d+$")
        & standardized["target"].str.match(r"^Gene_\d+$")
        & (standardized["source"] != standardized["target"])
    ]

    return standardized


def build_intermediate_networks(
    adata: "ad.AnnData",
    intermediate_dir: Path,
) -> None:
    """Build and save independently subsampled SCING edge lists."""
    number_of_samples, number_of_features = adata.shape

    neighbors_used = min(NNEIGHBORS, number_of_features - 1)
    pcs_used = min(NPCS, number_of_features - 1, number_of_samples - 1)

    if neighbors_used < 1 or pcs_used < 1:
        raise ValueError("The input matrix is too small for SCING.")

    intermediate_dir.mkdir(parents=True, exist_ok=True)

    for network_index in range(N_NETWORKS):
        grn = build.grnBuilder(
            adata=adata,
            ngenes=number_of_features,
            nneighbors=neighbors_used,
            npcs=pcs_used,
            subsample_perc=SUBSAMPLE_PERC,
            prefix=f"net.{network_index}",
            outdir=str(intermediate_dir),
            ncore=NCORE,
            mem_per_core=MEM_PER_CORE,
            verbose=False,
        )

        grn.subsample_cells()

        # The input is already feature-standardized. The additional
        # filter_genes transformation is omitted before connectivity filtering.
        intermediate_matrix = grn.adata.X
        if hasattr(intermediate_matrix, "toarray"):
            intermediate_matrix = intermediate_matrix.toarray()

        intermediate_matrix = np.asarray(intermediate_matrix, dtype=float)
        if not np.isfinite(intermediate_matrix).all():
            intermediate_matrix = np.nan_to_num(
                intermediate_matrix,
                nan=0.0,
                posinf=0.0,
                neginf=0.0,
            )

        grn.adata.X = intermediate_matrix
        grn.filter_gene_connectivities()
        grn.build_grn()
        grn.save_edges()


def read_intermediate_networks(
    intermediate_dir: Path,
) -> pd.DataFrame:
    """Read and standardize all intermediate SCING edge lists."""
    edge_tables = []

    for network_index in range(N_NETWORKS):
        filename = intermediate_dir / f"net.{network_index}.csv.gz"
        if not filename.is_file():
            raise FileNotFoundError(
                f"Intermediate SCING network not found: {filename}"
            )

        table = standardize_edge_columns(pd.read_csv(filename))
        table = (
            table.groupby(["source", "target"], as_index=False)
            .agg(importance=("importance", "mean"))
        )
        table["network_id"] = network_index
        edge_tables.append(table)

    return pd.concat(edge_tables, ignore_index=True)


# =============================================================================
# Consensus matrices
# =============================================================================

def build_consensus_edges(edge_table: pd.DataFrame) -> pd.DataFrame:
    """Aggregate directed edge weights and filter by appearance frequency."""
    consensus = (
        edge_table.groupby(["source", "target"], as_index=False)
        .agg(
            mean_importance=("importance", "mean"),
            n_appearance=("network_id", "nunique"),
        )
    )

    consensus["appearance_ratio"] = (
        consensus["n_appearance"] / float(N_NETWORKS)
    )

    return consensus[
        consensus["appearance_ratio"]
        >= MINIMUM_EDGE_APPEARANCE_THRESHOLD
    ].copy()


def directed_weight_matrix(
    consensus: pd.DataFrame,
    feature_names: np.ndarray,
) -> np.ndarray:
    """Convert directed consensus edges to a dense weight matrix."""
    feature_to_index = {
        feature: index for index, feature in enumerate(feature_names)
    }
    weights = np.zeros(
        (len(feature_names), len(feature_names)),
        dtype=float,
    )

    for row in consensus.itertuples(index=False):
        source_index = feature_to_index.get(str(row.source))
        target_index = feature_to_index.get(str(row.target))

        if (
            source_index is None
            or target_index is None
            or source_index == target_index
        ):
            continue

        weights[source_index, target_index] = float(row.mean_importance)

    np.fill_diagonal(weights, 0.0)
    return weights


def undirected_weight_matrix(
    directed_weights: np.ndarray,
    method: str,
) -> np.ndarray:
    """Combine opposite directed weights into a symmetric matrix."""
    if method == "max":
        undirected = np.maximum(directed_weights, directed_weights.T)
    elif method == "mean":
        undirected = (directed_weights + directed_weights.T) / 2.0
    elif method == "sum":
        undirected = directed_weights + directed_weights.T
    else:
        raise ValueError(
            "UNDIRECTED_COMBINE_METHOD must be 'max', 'mean', or 'sum'."
        )

    np.fill_diagonal(undirected, 0.0)
    return undirected


# =============================================================================
# Main analysis
# =============================================================================

def main() -> None:
    """Run SCING and save the directed and undirected consensus matrices."""
    warnings.filterwarnings("ignore")
    np.random.seed(RANDOM_SEED)

    os.environ["MKL_NUM_THREADS"] = str(NCORE)
    os.environ["NUMEXPR_NUM_THREADS"] = str(NCORE)
    os.environ["OMP_NUM_THREADS"] = str(NCORE)

    if N_NETWORKS < 1:
        raise ValueError("N_NETWORKS must be positive.")

    if not 0 < SUBSAMPLE_PERC <= 1:
        raise ValueError("SUBSAMPLE_PERC must be in (0, 1].")

    if not 0 <= MINIMUM_EDGE_APPEARANCE_THRESHOLD <= 1:
        raise ValueError(
            "MINIMUM_EDGE_APPEARANCE_THRESHOLD must be between 0 and 1."
        )

    # Initialize the directory containing the final SCING matrices.
    if OUTPUT_DIR.exists():
        shutil.rmtree(OUTPUT_DIR)
    OUTPUT_DIR.mkdir(parents=True)

    feature_by_sample = read_numeric_matrix(INPUT_MATRIX_FILE)
    feature_by_sample = subsample_samples(
        feature_by_sample,
        MAX_CELLS_FOR_SCING,
        RANDOM_SEED,
    )
    validate_feature_variance(feature_by_sample)

    number_of_features, number_of_samples = feature_by_sample.shape
    feature_names = generate_feature_names(number_of_features)
    sample_names = np.asarray(
        [f"Cell_{index + 1}" for index in range(number_of_samples)],
        dtype=str,
    )

    adata = ad.AnnData(X=feature_by_sample.T)
    adata.obs_names = sample_names
    adata.var_names = feature_names

    # Intermediate edge lists are temporary computational files.
    with tempfile.TemporaryDirectory(
        prefix="scing_intermediate_",
        dir=SCRIPT_DIR,
    ) as temporary_directory:
        intermediate_dir = Path(temporary_directory)

        build_intermediate_networks(
            adata=adata,
            intermediate_dir=intermediate_dir,
        )

        all_edges = read_intermediate_networks(intermediate_dir)
        consensus = build_consensus_edges(all_edges)

    directed = directed_weight_matrix(consensus, feature_names)
    undirected = undirected_weight_matrix(
        directed,
        UNDIRECTED_COMBINE_METHOD,
    )

    np.savetxt(
        DIRECTED_WEIGHT_MATRIX_FILE,
        directed,
        fmt="%.10g",
        delimiter="\t",
    )
    np.savetxt(
        UNDIRECTED_WEIGHT_MATRIX_FILE,
        undirected,
        fmt="%.10g",
        delimiter="\t",
    )


if __name__ == "__main__":
    main()
