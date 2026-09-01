#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
RegDiffusion inference for simulated local-feature expression matrices.

Each input matrix is interpreted as features by samples. The first N_LOCAL
rows are transposed to the samples-by-features format required by
RegDiffusion. The directed adjacency matrix and a symmetric local-local score
matrix are saved for each dataset.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np


# =============================================================================
# Analysis settings
# =============================================================================

# All relative paths are resolved from the directory containing this script.
SCRIPT_DIR = Path(__file__).resolve().parent

# Input matrices processed in sequence.
DATA_FILES = [
    "Data_whole_rep1.txt",
    "Data_whole_rep2.txt",
    "Data_whole_rep3.txt",
]

# Output directory corresponding to each input matrix.
OUTPUT_FOLDERS = [
    "RegDiffusion_Output_rep1",
    "RegDiffusion_Output_rep2",
    "RegDiffusion_Output_rep3",
]

# Number of local features retained from the beginning of each input matrix.
N_LOCAL = 500

# Number of RegDiffusion training steps.
STEPS = 1000

# Rule used to combine opposite directed scores.
# Supported values are "absmax", "absmean", and "signed_mean".
SYM_MODE = "absmax"

# Final output filenames.
DIRECTED_OUTPUT_NAME = "RegDiffusion_adj_directed.txt"
SYMMETRIC_OUTPUT_NAME = "RegDiffusion_LL_score_symmetric.txt"


# =============================================================================
# Matrix input and output
# =============================================================================

def read_numeric_matrix(filename: Path) -> np.ndarray:
    """Read a two-dimensional numeric matrix from TXT, TSV, CSV, or NPY."""
    if not filename.is_file():
        raise FileNotFoundError(f"Input matrix not found: {filename}")

    extension = filename.suffix.lower()

    if extension == ".npy":
        matrix = np.load(filename)
    elif extension == ".csv":
        matrix = np.loadtxt(filename, delimiter=",")
    elif extension == ".tsv":
        matrix = np.loadtxt(filename, delimiter="\t")
    elif extension in {".txt", ""}:
        matrix = np.loadtxt(filename)
    else:
        raise ValueError(
            "Supported input formats are TXT, TSV, CSV, and NPY."
        )

    matrix = np.asarray(matrix, dtype=np.float32)

    if matrix.ndim != 2:
        raise ValueError(
            f"Input matrix must be two-dimensional; received {matrix.shape}."
        )

    if not np.isfinite(matrix).all():
        raise ValueError("Input matrix contains NaN or infinite values.")

    return matrix


def write_matrix(filename: Path, matrix: np.ndarray) -> None:
    """Write a matrix as a tab-delimited text file."""
    np.savetxt(filename, matrix, fmt="%.16e", delimiter="\t")


# =============================================================================
# RegDiffusion input preparation
# =============================================================================

def prepare_regdiffusion_input(
    feature_by_sample: np.ndarray,
    number_of_local_features: int,
) -> np.ndarray:
    """Select local features and transpose the matrix to samples by features."""
    if number_of_local_features < 1:
        raise ValueError("N_LOCAL must be positive.")

    if feature_by_sample.shape[0] < number_of_local_features:
        raise ValueError(
            f"Input contains {feature_by_sample.shape[0]} feature rows, "
            f"but N_LOCAL is {number_of_local_features}."
        )

    local_feature_by_sample = feature_by_sample[
        :number_of_local_features,
        :,
    ]

    return local_feature_by_sample.T.copy()


# =============================================================================
# RegDiffusion inference
# =============================================================================

def run_regdiffusion_cpu(
    sample_by_feature: np.ndarray,
    training_steps: int,
) -> np.ndarray:
    """Fit RegDiffusion on CPU and return its directed adjacency matrix."""
    if training_steps < 1:
        raise ValueError("STEPS must be positive.")

    try:
        import regdiffusion as rd
    except ImportError as exc:
        raise ImportError(
            "The regdiffusion package must be installed before running "
            "this script."
        ) from exc

    try:
        trainer = rd.RegDiffusionTrainer(
            sample_by_feature,
            device="cpu",
        )
    except TypeError:
        trainer = rd.RegDiffusionTrainer(sample_by_feature)

    try:
        trainer.train(n_steps=training_steps)
    except TypeError:
        trainer.train()

    adjacency = np.asarray(trainer.get_adj(), dtype=np.float64)

    if (
        adjacency.ndim != 2
        or adjacency.shape[0] != adjacency.shape[1]
    ):
        raise ValueError(
            "RegDiffusion returned a non-square adjacency matrix: "
            f"{adjacency.shape}"
        )

    if not np.isfinite(adjacency).all():
        raise ValueError(
            "RegDiffusion returned NaN or infinite adjacency values."
        )

    np.fill_diagonal(adjacency, 0.0)
    return adjacency


def symmetric_local_score(
    directed_adjacency: np.ndarray,
    mode: str,
) -> np.ndarray:
    """Convert directed scores to a symmetric local-local score matrix."""
    if mode == "absmax":
        score = np.maximum(
            np.abs(directed_adjacency),
            np.abs(directed_adjacency.T),
        )
    elif mode == "absmean":
        score = 0.5 * (
            np.abs(directed_adjacency)
            + np.abs(directed_adjacency.T)
        )
    elif mode == "signed_mean":
        score = 0.5 * (
            directed_adjacency
            + directed_adjacency.T
        )
    else:
        raise ValueError(
            "SYM_MODE must be 'absmax', 'absmean', or 'signed_mean'."
        )

    np.fill_diagonal(score, 0.0)
    return score


# =============================================================================
# Dataset processing
# =============================================================================

def run_dataset(
    input_filename: str,
    output_folder: str,
) -> None:
    """Run RegDiffusion for one dataset and save the two final matrices."""
    input_path = SCRIPT_DIR / input_filename
    output_dir = SCRIPT_DIR / output_folder

    # Initialize the directory containing the final RegDiffusion matrices.
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True)

    feature_by_sample = read_numeric_matrix(input_path)
    sample_by_feature = prepare_regdiffusion_input(
        feature_by_sample,
        N_LOCAL,
    )

    directed_adjacency = run_regdiffusion_cpu(
        sample_by_feature,
        STEPS,
    )
    symmetric_score = symmetric_local_score(
        directed_adjacency,
        SYM_MODE,
    )

    write_matrix(
        output_dir / DIRECTED_OUTPUT_NAME,
        directed_adjacency,
    )
    write_matrix(
        output_dir / SYMMETRIC_OUTPUT_NAME,
        symmetric_score,
    )


# =============================================================================
# Main analysis
# =============================================================================

def main() -> None:
    """Process all configured datasets."""
    if len(DATA_FILES) != len(OUTPUT_FOLDERS):
        raise ValueError(
            "DATA_FILES and OUTPUT_FOLDERS must have the same length."
        )

    for input_filename, output_folder in zip(
        DATA_FILES,
        OUTPUT_FOLDERS,
    ):
        run_dataset(input_filename, output_folder)


if __name__ == "__main__":
    main()
