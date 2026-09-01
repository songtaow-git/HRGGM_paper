#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
GRNBoost2 inference for a numeric expression matrix.

The input matrix may be arranged as features by samples or samples by
features. GRNBoost2 estimates directed regulator-to-target importance scores,
which are converted to a weighted adjacency matrix for downstream analysis.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pandas as pd
from arboreto.algo import grnboost2
from distributed import Client, LocalCluster


# =============================================================================
# Analysis settings
# =============================================================================

# All relative paths are resolved from the directory containing this script.
SCRIPT_DIR = Path(__file__).resolve().parent

# Numeric input matrix without row or column labels.
INPUT_FILE = SCRIPT_DIR / "Data_whole_rep1.txt"

# Input orientation: "genes_by_cells" or "cells_by_genes".
MATRIX_ORIENTATION = "genes_by_cells"

# Directory containing the final weighted adjacency matrix.
OUTPUT_FOLDER = SCRIPT_DIR / "GRNBoost2_Output"

# Final output filename.
OUTPUT_FILENAME = "GRNBoost2_weighted_matrix.txt"

# Candidate-regulator configuration.
# When true, every feature is treated as a candidate regulator.
USE_ALL_GENES_AS_TF = True

# Optional list of generated feature labels, such as Gene_1 and Gene_5.
TF_LIST_FILE = SCRIPT_DIR / "TF_names.txt"

# When true, opposite directed weights are combined by their maximum value.
SYMMETRIZE_OUTPUT = True

# Dask resources used by GRNBoost2.
N_WORKERS = 4
THREADS_PER_WORKER = 1
USE_PROCESSES = False

# Random seed passed to GRNBoost2.
RANDOM_SEED = 123


# =============================================================================
# Input preparation
# =============================================================================

def read_numeric_matrix(filename: Path) -> np.ndarray:
    """Read a finite two-dimensional numeric matrix."""
    if not filename.is_file():
        raise FileNotFoundError(f"Input matrix not found: {filename}")

    matrix = np.loadtxt(filename, dtype=float)

    if matrix.ndim != 2:
        raise ValueError(
            f"Input matrix must be two-dimensional; received {matrix.shape}."
        )

    if not np.isfinite(matrix).all():
        raise ValueError("Input matrix contains NaN or infinite values.")

    return matrix


def prepare_expression_data(
    matrix: np.ndarray,
    orientation: str,
) -> tuple[pd.DataFrame, list[str]]:
    """Convert the input to the samples-by-features format used by GRNBoost2."""
    if orientation == "genes_by_cells":
        expression = matrix.T
        number_of_genes = matrix.shape[0]
        number_of_cells = matrix.shape[1]
    elif orientation == "cells_by_genes":
        expression = matrix
        number_of_cells = matrix.shape[0]
        number_of_genes = matrix.shape[1]
    else:
        raise ValueError(
            "MATRIX_ORIENTATION must be 'genes_by_cells' or "
            "'cells_by_genes'."
        )

    gene_names = [
        f"Gene_{index + 1}" for index in range(number_of_genes)
    ]
    cell_names = [
        f"Cell_{index + 1}" for index in range(number_of_cells)
    ]

    expression_data = pd.DataFrame(
        expression,
        index=cell_names,
        columns=gene_names,
    )

    return expression_data, gene_names


def load_candidate_regulators(
    gene_names: list[str],
) -> list[str]:
    """Return all features or a validated subset as candidate regulators."""
    if USE_ALL_GENES_AS_TF:
        return gene_names

    if not TF_LIST_FILE.is_file():
        raise FileNotFoundError(
            f"Candidate-regulator file not found: {TF_LIST_FILE}"
        )

    with TF_LIST_FILE.open("r", encoding="utf-8") as handle:
        requested_names = [
            line.strip() for line in handle if line.strip()
        ]

    valid_names = [
        name for name in requested_names if name in set(gene_names)
    ]

    if not valid_names:
        raise ValueError(
            "TF_LIST_FILE does not contain feature labels present in the "
            "input matrix."
        )

    return valid_names


# =============================================================================
# GRNBoost2 inference
# =============================================================================

def infer_grnboost2_network(
    expression_data: pd.DataFrame,
    regulator_names: list[str],
) -> pd.DataFrame:
    """Run GRNBoost2 and return directed regulator-to-target scores."""
    cluster = LocalCluster(
        n_workers=N_WORKERS,
        threads_per_worker=THREADS_PER_WORKER,
        processes=USE_PROCESSES,
        dashboard_address=None,
    )
    client = Client(cluster)

    try:
        network = grnboost2(
            expression_data=expression_data,
            tf_names=regulator_names,
            client_or_address=client,
            seed=RANDOM_SEED,
            verbose=False,
        )
    finally:
        client.close()
        cluster.close()

    required_columns = {"TF", "target", "importance"}
    if not required_columns.issubset(network.columns):
        raise ValueError(
            "GRNBoost2 returned an unexpected table structure: "
            f"{network.columns.tolist()}"
        )

    network = network.loc[
        network["TF"] != network["target"],
        ["TF", "target", "importance"],
    ].copy()

    network["importance"] = pd.to_numeric(
        network["importance"],
        errors="coerce",
    )
    network = network.dropna(subset=["importance"])

    return network


def build_weighted_matrix(
    network: pd.DataFrame,
    gene_names: list[str],
    symmetrize: bool,
) -> np.ndarray:
    """Convert directed importance scores to a dense weighted matrix."""
    gene_to_index = {
        gene: index for index, gene in enumerate(gene_names)
    }
    weighted_matrix = np.zeros(
        (len(gene_names), len(gene_names)),
        dtype=float,
    )

    for row in network.itertuples(index=False):
        regulator_index = gene_to_index.get(str(row.TF))
        target_index = gene_to_index.get(str(row.target))

        if (
            regulator_index is None
            or target_index is None
            or regulator_index == target_index
        ):
            continue

        importance = float(row.importance)
        if importance > weighted_matrix[regulator_index, target_index]:
            weighted_matrix[regulator_index, target_index] = importance

    if symmetrize:
        weighted_matrix = np.maximum(
            weighted_matrix,
            weighted_matrix.T,
        )

    np.fill_diagonal(weighted_matrix, 0.0)
    return weighted_matrix


# =============================================================================
# Main analysis
# =============================================================================

def main() -> None:
    """Run GRNBoost2 and save the weighted adjacency matrix."""
    if N_WORKERS < 1 or THREADS_PER_WORKER < 1:
        raise ValueError(
            "N_WORKERS and THREADS_PER_WORKER must be positive."
        )

    # Initialize the directory containing the final GRNBoost2 matrix.
    if OUTPUT_FOLDER.exists():
        shutil.rmtree(OUTPUT_FOLDER)
    OUTPUT_FOLDER.mkdir(parents=True)

    input_matrix = read_numeric_matrix(INPUT_FILE)
    expression_data, gene_names = prepare_expression_data(
        input_matrix,
        MATRIX_ORIENTATION,
    )
    regulator_names = load_candidate_regulators(gene_names)

    network = infer_grnboost2_network(
        expression_data,
        regulator_names,
    )
    weighted_matrix = build_weighted_matrix(
        network,
        gene_names,
        SYMMETRIZE_OUTPUT,
    )

    np.savetxt(
        OUTPUT_FOLDER / OUTPUT_FILENAME,
        weighted_matrix,
        fmt="%.10g",
        delimiter="\t",
    )


if __name__ == "__main__":
    main()
