################################################################################
## Thresholded Graphical Lasso (TGL)
## Wang and Allen, Biometrika 2023
##
## The script fits graphical lasso at a small regularization level and applies
## hard thresholding to the off-diagonal precision-matrix entries. The threshold
## is selected to retain a specified number of undirected local-local edges.
##
## Outputs:
##   TGL_Theta_LL_initial_small_lambda_glasso.txt
##       Initial graphical-lasso precision matrix before hard thresholding.
##
##   TGL_Theta_LL_primary_thresholded.txt
##       Final thresholded local-local precision matrix.
################################################################################

rm(list = ls())

options(stringsAsFactors = FALSE)
options(scipen = 999)

## Directory containing the two final precision matrices.
OUTPUT_DIR <- file.path(SCRIPT_DIR, "TGL_Output_rep1")

## Pure-numeric input matrix with no row or column names.
INPUT_FILE <- file.path(SCRIPT_DIR, "Data_whole_rep1.txt")

################################################################################
## Analysis settings
################################################################################

get_script_directory <- function() {
  command_arguments <- commandArgs(trailingOnly = FALSE)
  file_argument <- grep("^--file=", command_arguments, value = TRUE)

  if (length(file_argument) > 0L) {
    script_path <- sub("^--file=", "", file_argument[1L])
    return(dirname(normalizePath(
      script_path,
      winslash = "/",
      mustWork = TRUE
    )))
  }

  frame_files <- vapply(
    sys.frames(),
    function(frame) {
      if (is.null(frame$ofile)) "" else as.character(frame$ofile)
    },
    FUN.VALUE = character(1)
  )
  frame_files <- frame_files[nzchar(frame_files)]

  if (length(frame_files) > 0L) {
    return(dirname(normalizePath(
      frame_files[length(frame_files)],
      winslash = "/",
      mustWork = TRUE
    )))
  }

  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

SCRIPT_DIR <- get_script_directory()

## Matrix orientation in INPUT_FILE:
##   "features_by_samples" for features x samples
##   "samples_by_features" for samples x features
INPUT_ORIENTATION <- "features_by_samples"

## Number of undirected local-local edges retained after hard thresholding.
TRUE_UNDIRECTED_EDGE_COUNT <- 600L

## Absolute tolerance used to count nonzero edges in the initial solution.
NUMERICAL_ZERO <- 1e-12

## Initial graphical-lasso penalty:
##   lambda0 = INITIAL_LAMBDA_MULTIPLIER * sqrt(log(p) / n)
INITIAL_LAMBDA_MULTIPLIER <- 1.0

## Maximum number of deterministic halvings applied when the initial solution
## contains fewer edges than TRUE_UNDIRECTED_EDGE_COUNT.
MAX_LAMBDA_HALVINGS <- 12L

################################################################################
## Packages
################################################################################

required_packages <- c("data.table", "glasso")
missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0L) {
  stop(
    "Required R packages are not installed: ",
    paste(missing_packages, collapse = ", ")
  )
}

################################################################################
## Input and output validation
################################################################################

if (!file.exists(INPUT_FILE)) {
  stop("Input file does not exist: ", INPUT_FILE)
}

if (!INPUT_ORIENTATION %in% c(
  "features_by_samples",
  "samples_by_features"
)) {
  stop(
    "INPUT_ORIENTATION must be 'features_by_samples' or ",
    "'samples_by_features'."
  )
}

if (!is.numeric(TRUE_UNDIRECTED_EDGE_COUNT) ||
    length(TRUE_UNDIRECTED_EDGE_COUNT) != 1L ||
    !is.finite(TRUE_UNDIRECTED_EDGE_COUNT) ||
    TRUE_UNDIRECTED_EDGE_COUNT < 1) {
  stop("TRUE_UNDIRECTED_EDGE_COUNT must be one positive integer.")
}

TRUE_UNDIRECTED_EDGE_COUNT <- as.integer(TRUE_UNDIRECTED_EDGE_COUNT)

if (!is.numeric(INITIAL_LAMBDA_MULTIPLIER) ||
    length(INITIAL_LAMBDA_MULTIPLIER) != 1L ||
    !is.finite(INITIAL_LAMBDA_MULTIPLIER) ||
    INITIAL_LAMBDA_MULTIPLIER <= 0) {
  stop("INITIAL_LAMBDA_MULTIPLIER must be one positive number.")
}

if (!is.numeric(MAX_LAMBDA_HALVINGS) ||
    length(MAX_LAMBDA_HALVINGS) != 1L ||
    !is.finite(MAX_LAMBDA_HALVINGS) ||
    MAX_LAMBDA_HALVINGS < 0) {
  stop("MAX_LAMBDA_HALVINGS must be one nonnegative integer.")
}

MAX_LAMBDA_HALVINGS <- as.integer(MAX_LAMBDA_HALVINGS)

## Initialize a clean output directory for the final TGL matrices.
if (dir.exists(OUTPUT_DIR)) {
  unlink(OUTPUT_DIR, recursive = TRUE, force = TRUE)
}

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

################################################################################
## Matrix utilities
################################################################################

count_undirected_edges <- function(theta, zero_tolerance = 0) {
  if (!is.matrix(theta) || nrow(theta) != ncol(theta)) {
    stop("theta must be a square matrix.")
  }

  sum(abs(theta[upper.tri(theta)]) > zero_tolerance)
}

hard_threshold_off_diagonal <- function(theta_initial, threshold) {
  theta_thresholded <- theta_initial
  off_diagonal <- row(theta_thresholded) != col(theta_thresholded)

  theta_thresholded[
    off_diagonal & abs(theta_thresholded) <= threshold
  ] <- 0

  (theta_thresholded + t(theta_thresholded)) / 2
}

oracle_threshold_from_edge_count <- function(
    theta_initial,
    target_edge_count
) {
  absolute_weights <- abs(theta_initial[upper.tri(theta_initial)])

  if (any(!is.finite(absolute_weights))) {
    stop(
      "The initial precision matrix contains non-finite ",
      "off-diagonal entries."
    )
  }

  total_possible_edges <- length(absolute_weights)

  if (target_edge_count > total_possible_edges) {
    stop(
      "TRUE_UNDIRECTED_EDGE_COUNT exceeds the number of possible ",
      "undirected edges."
    )
  }

  sorted_weights <- sort(absolute_weights, decreasing = TRUE)

  if (target_edge_count == total_possible_edges) {
    return(-Inf)
  }

  kth_weight <- sorted_weights[target_edge_count]
  next_weight <- sorted_weights[target_edge_count + 1L]

  if (kth_weight <= 0) {
    stop(
      "The initial graphical-lasso solution contains fewer than ",
      target_edge_count,
      " nonzero undirected edges."
    )
  }

  if (kth_weight == next_weight) {
    stop(
      "A tie occurs at the edge-count boundary; a scalar threshold cannot ",
      "retain exactly ",
      target_edge_count,
      " edges."
    )
  }

  (kth_weight + next_weight) / 2
}

write_matrix_tab <- function(matrix, filename) {
  write.table(
    matrix,
    file = filename,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
}

################################################################################
## Read and standardize the input matrix
################################################################################

input_table <- data.table::fread(
  INPUT_FILE,
  header = FALSE,
  data.table = FALSE,
  showProgress = FALSE
)

input_matrix <- as.matrix(input_table)
storage.mode(input_matrix) <- "double"

if (any(!is.finite(input_matrix))) {
  stop("Input data contain nonnumeric or non-finite values.")
}

if (INPUT_ORIENTATION == "features_by_samples") {
  X <- t(input_matrix)
} else {
  X <- input_matrix
}

rm(input_table, input_matrix)

n <- nrow(X)
p <- ncol(X)

if (n < 2L || p < 2L) {
  stop("The analysis matrix must contain at least two samples and features.")
}

maximum_possible_edges <- p * (p - 1L) / 2L
if (TRUE_UNDIRECTED_EDGE_COUNT > maximum_possible_edges) {
  stop(
    "TRUE_UNDIRECTED_EDGE_COUNT exceeds p * (p - 1) / 2."
  )
}

feature_means <- colMeans(X)
feature_sds <- apply(X, 2, sd)

if (any(!is.finite(feature_sds)) || any(feature_sds <= 0)) {
  stop("The input contains constant or invalid features.")
}

X <- sweep(X, 2, feature_means, FUN = "-")
X <- sweep(X, 2, feature_sds, FUN = "/")

sample_correlation <- crossprod(X) / n
sample_correlation <- (sample_correlation + t(sample_correlation)) / 2
diag(sample_correlation) <- 1

rm(X, feature_means, feature_sds)

################################################################################
## Initial small-lambda graphical lasso
################################################################################

paper_rate <- sqrt(log(p) / n)
lambda_multiplier <- INITIAL_LAMBDA_MULTIPLIER
lambda0 <- lambda_multiplier * paper_rate

theta_initial <- NULL
initial_edge_count <- NA_integer_

for (halving_step in 0:MAX_LAMBDA_HALVINGS) {
  fit_initial <- tryCatch(
    glasso::glasso(
      s = sample_correlation,
      rho = lambda0,
      penalize.diagonal = FALSE
    ),
    error = function(e) NULL
  )

  if (!is.null(fit_initial)) {
    theta_initial <- (fit_initial$wi + t(fit_initial$wi)) / 2

    initial_edge_count <- count_undirected_edges(
      theta_initial,
      zero_tolerance = NUMERICAL_ZERO
    )

    if (initial_edge_count >= TRUE_UNDIRECTED_EDGE_COUNT) {
      break
    }
  }

  lambda_multiplier <- lambda_multiplier / 2
  lambda0 <- lambda_multiplier * paper_rate
}

if (is.null(theta_initial)) {
  stop(
    "The graphical-lasso fit failed for every attempted initial penalty."
  )
}

if (initial_edge_count < TRUE_UNDIRECTED_EDGE_COUNT) {
  stop(
    "The densest attempted graphical-lasso solution contains fewer than ",
    TRUE_UNDIRECTED_EDGE_COUNT,
    " undirected edges."
  )
}

################################################################################
## Hard thresholding to the specified edge count
################################################################################

tau <- oracle_threshold_from_edge_count(
  theta_initial = theta_initial,
  target_edge_count = TRUE_UNDIRECTED_EDGE_COUNT
)

theta_tgl <- hard_threshold_off_diagonal(
  theta_initial = theta_initial,
  threshold = tau
)

final_edge_count <- count_undirected_edges(
  theta_tgl,
  zero_tolerance = 0
)

if (final_edge_count != TRUE_UNDIRECTED_EDGE_COUNT) {
  stop(
    "The thresholded matrix contains ",
    final_edge_count,
    " edges instead of ",
    TRUE_UNDIRECTED_EDGE_COUNT,
    "."
  )
}

################################################################################
## Save final matrices
################################################################################

write_matrix_tab(
  theta_initial,
  file.path(
    OUTPUT_DIR,
    "TGL_Theta_LL_initial_small_lambda_glasso.txt"
  )
)

write_matrix_tab(
  theta_tgl,
  file.path(
    OUTPUT_DIR,
    "TGL_Theta_LL_primary_thresholded.txt"
  )
)
