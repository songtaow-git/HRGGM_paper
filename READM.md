# HR-GGM paper analysis and reproducibility code

This repository contains the code used to generate and evaluate the results reported in the HR-GGM manuscript. The analyses cover:

- recovery of the network among observed variables in simulation;
- recovery and specificity of associations between the target hidden master regulator H1 and observed variables;
- comparison of Ca$^{2+}$-associated gene rankings in mouse astrocyte single-cell RNA-seq data;
- supplementary analyses of noise sensitivity, observed-network sparsity, and hyperparameter selection.

The code archive and the data archive are distributed separately:

- `HRGGM_paper.zip`: source code, analysis scripts, and representative output files;
- `Data_HRGGM_paper.zip`: input data, fitted matrices, and intermediate results required by the analysis scripts.

Download `Data_HRGGM_paper.zip` from the data repository associated with the paper: **[insert the Zenodo record or DOI URL here]**.

## 1. Required directory assembly

After extraction, the two archives contain the same five top-level analysis directories:

```text
HRGGM_paper/
├── Paper_astrocyte_data_method_comparison/
├── Paper_simulation_data_method_comparison/
├── Supplementary_simulation_analysis_S1/
├── Supplementary_simulation_analysis_S2/
└── Supplementary_simulation_analysis_S3/

Data_HRGGM_paper/
├── Paper_astrocyte_data_method_comparison/
├── Paper_simulation_data_method_comparison/
├── Supplementary_simulation_analysis_S1/
├── Supplementary_simulation_analysis_S2/
└── Supplementary_simulation_analysis_S3/
```

Copy the **contents** of each directory under `Data_HRGGM_paper` into the identically named directory under `HRGGM_paper`. Preserve every internal directory name and relative path.

The final working copy must have one combined tree rooted at `HRGGM_paper`. Do not place the entire data archive as an additional nested directory, and do not flatten individual data files into the repository root.

### Windows File Explorer

1. Extract both ZIP archives into the same parent directory.
2. Open `Data_HRGGM_paper` and copy its five top-level directories.
3. Paste them into `HRGGM_paper`.
4. Allow Windows to merge directories with identical names.

The equivalent PowerShell command, run from the common parent directory, is:

```powershell
robocopy .\Data_HRGGM_paper .\HRGGM_paper /E
if ($LASTEXITCODE -ge 8) { throw "Data copy failed with robocopy exit code $LASTEXITCODE" }
```

Do not add the `/MIR` option, because it can delete code files that are absent from the data archive.

### macOS or Linux

From the common parent directory, run:

```bash
rsync -a Data_HRGGM_paper/ HRGGM_paper/
```

After merging, a module should contain its original scripts together with the data and fitted-result directories expected by those scripts. For example:

```text
HRGGM_paper/Paper_simulation_data_method_comparison/
├── C_code_HRGGM/
├── Code_published_methods/
├── Data/
├── GT/
├── HRGGM/
├── LVGGM/
├── Prior_genes/
├── Code_network_among_observed_variables.m
└── Code_edges_between_H1_and_observed_variables.m
```

## 2. Terminology

The manuscript uses observed variables (`O`), hidden master regulatory variables (`H`), and hidden global variables (`G`). Some internal variable names and simulation directory identifiers retain earlier notation for compatibility:

| Legacy code notation | Current manuscript terminology |
|---|---|
| `local`, `L` | observed variable, `O` |
| `LL` | observed--observed block or edge, `O--O` |
| `LH` | observed--hidden-regulatory block or edge, `O--H` |
| `LG` | observed--global block or edge, `O--G` |

Do not rename legacy input fields or dataset directories. Several scripts intentionally recognize names such as `LLonly`, `LLbetween`, `Actual_LL_total_edges`, and `Data_200L_LL250_LH040_*` while reporting outputs with the current terminology.

## 3. Software requirements

The primary paper figures and tables are generated in MATLAB. The following environment is recommended:

- MATLAB R2020a or later;
- Statistics and Machine Learning Toolbox, used by functions including `corr`, `fishertest`, and `chi2cdf`;
- a spreadsheet writer supported by MATLAB for the `.xlsx` outputs.

Optional recomputation of model and comparison-method outputs additionally requires:

- a C++17 compiler, CMake 3.16 or later, and OpenMP for HR-GGM; Eigen headers are bundled under each `C_code_HRGGM` directory;
- R with `data.table` and `glasso` for TGL;
- Python 3 with the dependencies required by GRNBoost2, SCING, and RegDiffusion;
- the bundled LogdetPPA code and compatible MEX binaries, or locally rebuilt MEX files, for LV-GGM.

Third-party method wrappers are provided for traceability. Their configuration blocks specify the input repeat, paths, resources, and output directory. Direct figure reproduction does not require rerunning these methods when the archived fitted matrices have been merged from `Data_HRGGM_paper.zip`.

## 4. Analysis modules

| Directory | Manuscript analysis | Primary script(s) | Principal output |
|---|---|---|---|
| `Paper_simulation_data_method_comparison` | Main-text simulation: network among observed variables | `Code_network_among_observed_variables.m` | `Observed_variable_network_analysis/Observed_variable_network_recovery.png` and `Observed_variable_network_metrics.xlsx` |
| `Paper_simulation_data_method_comparison` | Main-text simulation: H1--observed associations and signal specificity | `Code_edges_between_H1_and_observed_variables.m` | `H1_observed_edge_analysis/H1_observed_edge_recovery_and_signal_specificity_1x3.png` and `H1_observed_edge_analysis.xlsx` |
| `Paper_astrocyte_data_method_comparison` | Main-text astrocyte comparison of Ca$^{2+}$-associated gene rankings | `Code_Ca_method_comparison.m`; `Code_statistic_test.m` | `Method_Comparison_Output/MethodComparison_Main_ABC_1x3.png`, ranking tables, resampling results, and threshold-based enrichment statistics |
| `Supplementary_simulation_analysis_S1` | Robustness to noise and fitting hidden variables when true hidden effects are absent | `Code_simulation_sensitivity.m` | `Sensitive_O_OH_Analysis/Figure_O_OH_combined_reordered_CleanLines.png` and threshold/error tables |
| `Supplementary_simulation_analysis_S2` | Recovery of observed-variable edges across network sparsity levels | `Code_HRGGM_O_sparsity.m` | `O_Sparsity_Performance_Output/Figure_O_sparsity_performance_1x3.png` and condition-level metrics |
| `Supplementary_simulation_analysis_S3` | Post hoc validation of data-driven selection of $\lambda$ and $\alpha$ | `Code_parameter_selection_validation.m` | `Parameter_Selection_Validation_OUTPUT/Figure_Main_ParameterSelection_1x3.png`, scan tables, and validation summaries |

### 4.1 Main-text simulation comparison

Set the MATLAB current folder to `Paper_simulation_data_method_comparison` before running the scripts.

`Code_network_among_observed_variables.m` evaluates all unordered pairs in the 500-by-500 observed-variable block across three matched repeats. It compares HR-GGM, LV-GGM, TGL, GRNBoost2, SCING, and RegDiffusion using absolute undirected edge scores, precision--recall curves, and candidate-list recovery.

`Code_edges_between_H1_and_observed_variables.m` evaluates 45 held-out H1 targets among 470 non-prior observed variables. It compares HR-GGM, a prior-selected LV-GGM component, and prior correlation. The same script also calculates the regression-based composition of the recovered H1 signal; the LV-GGM component used for this composition analysis is selected by an oracle ground-truth criterion, as specified in the manuscript.

The optional `Data_generation.m` script recreates the simulated data design with 500 observed variables, five ground-truth hidden regulatory variables, one hidden global variable, and three repeats. The archived data should be used for exact reproduction of the reported results.

### 4.2 Astrocyte comparison

Set the MATLAB current folder to `Paper_astrocyte_data_method_comparison`.

`Code_Ca_method_comparison.m` requires the following files in that directory:

```text
S3_Precision_matrix.txt
S1_whole_gene_list.xlsx
S5_Ca_associated_reference_label.xlsx
LVGGM_L_lowrank.txt
Data_whole.txt
```

The script ranks the 2,958 non-prior genes using the HR-GGM observed--hidden-regulatory precision entries, a prior-selected LV-GGM component, and mean absolute correlation with the 42 prior genes. It produces top-$K$ positive rates, enrichment folds, precision--recall results, permutation tests, bootstrap comparisons, and spectral diagnostics.

`Code_statistic_test.m` reproduces the threshold-based contingency analysis at an absolute HR-GGM score of 0.038 and reports the odds ratio, 95% confidence interval, Fisher's exact test, and chi-square test in the MATLAB command window.

### 4.3 Supplementary analysis S1

Set the MATLAB current folder to `Supplementary_simulation_analysis_S1` and run `Code_simulation_sensitivity.m`. The script discovers the archived condition directories beside itself, reads matched `Theta_full.txt` and `Theta_final.txt` matrices, and evaluates:

- observed-network support recovery with and without true hidden effects;
- the threshold range retaining at least 90% of the maximum F1 score;
- relative Frobenius error of the observed block;
- observed--hidden-regulatory recovery and sign-invariant coefficient error when true regulatory effects are present.

Legacy condition names containing `LLbetween` or `LLonly` must remain unchanged because they are used during dataset discovery.

### 4.4 Supplementary analysis S2

Set the MATLAB current folder to `Supplementary_simulation_analysis_S2` and run `Code_HRGGM_O_sparsity.m`. The script reads `All_conditions_summary.csv` and the corresponding condition directories. For each condition it expects the ground truth and estimate at paths equivalent to:

```text
<condition>/build/Data/Theta_full.txt
<condition>/build/Result/Theta_final.txt
```

It summarizes normalized AUPRC, tie-aware top-$K$ recovery, and relative off-diagonal Frobenius error over 36 datasets spanning four observed-network densities, three noise levels, and three repeats.

### 4.5 Supplementary analysis S3

Set the MATLAB current folder to `Supplementary_simulation_analysis_S3` and run `Code_parameter_selection_validation.m`. The script discovers directories matching the archived pattern:

```text
Data_200L_LL250_LH040_noise_*pct_rep_*
```

Within each dataset, it reads the candidate grids, the selected-score histories, the fitted matrices, the full-data parameter scans, and the ground-truth matrix. Ground truth is used only for post hoc validation. The output compares the selected parameters with favorable ground-truth-based reference regions for observed-block error and the population analogue of the model-fit/adjusted-entropy score.

## 5. Recommended execution order

For direct reproduction from the archived inputs and fitted matrices:

1. Merge `Data_HRGGM_paper` into `HRGGM_paper` as described in Section 1.
2. Start MATLAB and set the current folder to one analysis module at a time.
3. Run the two main-text simulation analysis scripts.
4. Run `Code_Ca_method_comparison.m`, followed by `Code_statistic_test.m`.
5. Run the S1, S2, and S3 analysis scripts independently.

These modules do not need to be run in a single MATLAB session. Each primary analysis script writes to its own output directory.

To reproduce the complete computational pipeline from generated observations rather than from archived fitted matrices, run the relevant `Data_generation.m`, fit HR-GGM with the module-specific `C_code_HRGGM`, run the required comparison methods, and then run the MATLAB analysis script. Read `C_code_HRGGM/README_Data_and_Main_mapping.txt` for the HR-GGM input-file format. Each module's `Main.cpp` is the authoritative run-specific configuration; do not replace it with a copy from another module.

## 6. Important reproducibility notes

- Keep all matrix row and column orders unchanged.
- C++ index files use zero-based indices; the MATLAB analysis scripts convert them to one-based indices where required.
- `Theta_full.txt` is simulation ground truth and is not used by HR-GGM parameter selection or fitting.
- The archived fitted matrices and intermediate score files are the inputs for exact figure reproduction.
- Several scripts overwrite files in their output directories. `Code_parameter_selection_validation.m` recreates its output directory, and the simulation generators have `reset_output_folder = true`. Preserve a copy of archived outputs before intentionally regenerating data or results.
- Do not rename source dataset directories or legacy table columns solely to match manuscript terminology; the analysis scripts already translate exported descriptions where appropriate.

## 7. Citation and third-party software

If you use this code or the accompanying data, please cite the HR-GGM paper:

> **HR-GGM paper citation to be added after publication.**

The comparison directory contains code or wrappers associated with LV-GGM, TGL, GRNBoost2, SCING, and RegDiffusion. Please also cite the corresponding original methods when using or rerunning them. The bundled Eigen and LogdetPPA components retain their own license and citation files where provided.

## 8. Data and code availability

- Code: this GitHub repository (`HRGGM_paper`)
- Data and archived fitted results: `Data_HRGGM_paper.zip` on Zenodo (**insert DOI URL**)

The code and data archives are versioned separately. For reproducibility, report both the GitHub commit or release and the Zenodo record version used in an analysis.
