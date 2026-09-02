# HR-GGM paper analysis and reproducibility code

This repository contains the code used to generate and evaluate the results reported in the HR-GGM manuscript. The analyses cover:

* recovery of the network among observed variables in simulation;
* recovery and specificity of associations between the target hidden master regulator H1 and observed variables;
* comparison of Ca²⁺-associated gene rankings in mouse astrocyte single-cell RNA-seq data;
* supplementary analyses of recovery robustness, observed-network sparsity, and parameter selection.

The analysis code is stored directly in this GitHub repository. The data and archived fitted results are distributed separately:

* `Data\\\_HRGGM\\\_paper.zip`: input data, fitted matrices, and intermediate results required by the analysis scripts.

Download [`Data\\\_HRGGM\\\_paper.zip` from the Zenodo draft record](https://zenodo.org/records/22217595?preview=1&token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjMwZjVkZTQ3LTAxMzYtNDEyNy04ODllLTUzZTRkNjQ5ZTU2MiIsImRhdGEiOnt9LCJyYW5kb20iOiI5NmU1MTNhMGQyNjA3ZTNlY2NhNzhkNTgzYzliNmQ0MSJ9.4h8NHWIfglk8OPAE7i7CWuAb8tDwtuTK2Bnncz2Y7SdJj_6rkwYNO_lcI6AqRKLJzzje_dmQGqy6YtR--EteSg). The current link provides preview access before the Zenodo record is formally published.

The repository also contains a `Supplementary\\\_tables/` directory. These files
are provided directly through GitHub and are not included in
`Data\\\_HRGGM\\\_paper.zip`:

* `S1\\\_whole\\\_gene\\\_list`: complete gene list used in the astrocyte analysis;
* `S2\\\_Network\\\_overview`: summary of the estimated HR-GGM network;
* `S3\\\_Precision\\\_matrix`: archived estimated precision matrix;
* `S4\\\_Enrichr\\\_result`: Enrichr functional-enrichment results;
* `S5\\\_Ca\\\_associated\\\_reference\\\_label`: reference labels for Ca²⁺-associated genes;
* `S6\\\_Ca2\\\_literature\\\_evidence`: literature evidence for Ca²⁺-associated genes.

## 1\. Required directory assembly

The GitHub repository and the extracted data archive contain the same five top-level analysis directories:

```text
HRGGM\\\_paper/
├── Paper\\\_astrocyte\\\_data\\\_method\\\_comparison/
├── Paper\\\_simulation\\\_data\\\_method\\\_comparison/
├── Supplementary\\\_simulation\\\_analysis\\\_S1/
├── Supplementary\\\_simulation\\\_analysis\\\_S2/
├── Supplementary\\\_simulation\\\_analysis\\\_S3/
└── Supplementary\\\_tables/                  # GitHub only

Data\\\_HRGGM\\\_paper/
├── Paper\\\_astrocyte\\\_data\\\_method\\\_comparison/
├── Paper\\\_simulation\\\_data\\\_method\\\_comparison/
├── Supplementary\\\_simulation\\\_analysis\\\_S1/
├── Supplementary\\\_simulation\\\_analysis\\\_S2/
└── Supplementary\\\_simulation\\\_analysis\\\_S3/
```

Copy the **contents** of each directory under `Data\\\_HRGGM\\\_paper` into the identically named directory under `HRGGM\\\_paper`. Preserve every internal directory name and relative path.

`Supplementary\\\_tables/` already belongs to the GitHub repository and has no
matching directory in the Zenodo archive. Leave it unchanged when merging the
five analysis directories.

The final working copy must have one combined tree rooted at `HRGGM\\\_paper`. Do not place the entire data archive as an additional nested directory, and do not flatten individual data files into the repository root.

### Windows File Explorer

1. Clone or download this GitHub repository as the folder `HRGGM\\\_paper`.
2. Extract `Data\\\_HRGGM\\\_paper.zip` beside the repository folder.
3. Open `Data\\\_HRGGM\\\_paper` and copy its five top-level directories.
4. Paste them into `HRGGM\\\_paper` and allow Windows to merge directories with identical names.

The equivalent PowerShell command, run from the common parent directory, is:

```powershell
robocopy .\\\\Data\\\_HRGGM\\\_paper .\\\\HRGGM\\\_paper /E
if ($LASTEXITCODE -ge 8) { throw "Data copy failed with robocopy exit code $LASTEXITCODE" }
```

Do not add the `/MIR` option, because it can delete code files that are absent from the data archive.

### macOS or Linux

From the common parent directory, run:

```bash
rsync -a Data\\\_HRGGM\\\_paper/ HRGGM\\\_paper/
```

After merging, each module contains its original code together with the corresponding data and archived fitted results. For example, the complete top level of the merged simulation-comparison module is:

```text
HRGGM\\\_paper/Paper\\\_simulation\\\_data\\\_method\\\_comparison/
├── C\\\_code\\\_HRGGM/
├── Code\\\_published\\\_methods/
├── Data/
├── GRNBoost2/
├── GT/
├── H1\\\_observed\\\_edge\\\_analysis/
├── HRGGM/
├── LVGGM/
├── Observed\\\_variable\\\_network\\\_analysis/
├── Prior\\\_genes/
├── RegDiffusion/
├── SCING/
├── TGL/
├── Alpha\\\_list.txt
├── Code\\\_edges\\\_between\\\_H1\\\_and\\\_observed\\\_variables.m
├── Code\\\_network\\\_among\\\_observed\\\_variables.m
├── Data\\\_generation.m
└── Lambda\\\_list.txt
```

## 2\. Terminology

The manuscript uses observed variables (`O`), hidden master regulatory variables (`H`), and hidden global variables (`G`). Some internal variable names and simulation directory identifiers retain earlier notation for compatibility:

|Legacy code notation|Current manuscript terminology|
|-|-|
|`local`, `L`|observed variable, `O`|
|`LL`|observed--observed block or edge, `O--O`|
|`LH`|observed--hidden-regulatory block or edge, `O--H`|
|`LG`|observed--global block or edge, `O--G`|

Do not rename legacy input fields or dataset directories. Several scripts intentionally recognize names such as `LLonly`, `LLbetween`, `Actual\\\_LL\\\_total\\\_edges`, and `Data\\\_200L\\\_LL250\\\_LH040\\\_\\\*` while reporting outputs with the current terminology.

## 3\. Software requirements

The primary paper figures and tables are generated in MATLAB. The following environment is recommended:

* MATLAB R2020a or later;
* Statistics and Machine Learning Toolbox, used by functions including `corr`, `fishertest`, and `chi2cdf`;
* a spreadsheet writer supported by MATLAB for the `.xlsx` outputs.

Optional recomputation of model and comparison-method outputs additionally requires:

* a C++17 compiler, CMake 3.16 or later, and OpenMP for HR-GGM; Eigen headers are bundled under each `C\\\_code\\\_HRGGM` directory;
* R with `data.table` and `glasso` for TGL;
* Python 3 with the dependencies required by GRNBoost2, SCING, and RegDiffusion;
* the bundled LogdetPPA code and compatible MEX binaries, or locally rebuilt MEX files, for LV-GGM.

Third-party method wrappers are provided for traceability. Their configuration blocks specify the input repeat, paths, resources, and output directory. Direct figure reproduction does not require rerunning these methods when the archived fitted matrices have been merged from `Data\\\_HRGGM\\\_paper.zip`.

## 4\. Analysis modules

The repository is organized by manuscript analysis:

* **Main-text simulation comparison**  
Directory: `Paper\\\_simulation\\\_data\\\_method\\\_comparison`  
Scripts: `Code\\\_network\\\_among\\\_observed\\\_variables.m` and `Code\\\_edges\\\_between\\\_H1\\\_and\\\_observed\\\_variables.m`
* **Main-text astrocyte comparison of Ca²⁺-associated gene rankings**  
Directory: `Paper\\\_astrocyte\\\_data\\\_method\\\_comparison`  
Scripts: `Code\\\_Ca\\\_method\\\_comparison.m` and `Code\\\_statistic\\\_test.m`
* **Supplementary S1: Robustness of observed-variable network and hidden-regulator association recovery**  
Directory: `Supplementary\\\_simulation\\\_analysis\\\_S1`  
Script: `Code\\\_simulation\\\_sensitivity.m`
* **Supplementary S2: Recovery of observed-variable edges under different sparsity levels**  
Directory: `Supplementary\\\_simulation\\\_analysis\\\_S2`  
Script: `Code\\\_HRGGM\\\_O\\\_sparsity.m`
* **Supplementary S3: Validation of the data-driven parameter-selection procedure**  
Directory: `Supplementary\\\_simulation\\\_analysis\\\_S3`  
Script: `Code\\\_parameter\\\_selection\\\_validation.m`

### 4.1 Paper simulation data experiment

Set the MATLAB current folder to `Paper\\\_simulation\\\_data\\\_method\\\_comparison` before running the scripts.

`Code\\\_network\\\_among\\\_observed\\\_variables.m` evaluates all unordered pairs in the 500-by-500 observed-variable block across three matched repeats. It compares HR-GGM, LV-GGM, TGL, GRNBoost2, SCING, and RegDiffusion using absolute undirected edge scores, precision--recall curves, and candidate-list recovery.

`Code\\\_edges\\\_between\\\_H1\\\_and\\\_observed\\\_variables.m` evaluates 45 held-out H1 targets among 470 non-prior observed variables. It compares HR-GGM, a prior-selected LV-GGM component, and prior correlation. The same script also calculates the regression-based composition of the recovered H1 signal; the LV-GGM component used for this composition analysis is selected by an oracle ground-truth criterion, as specified in the manuscript.

The principal outputs are `Observed\\\_variable\\\_network\\\_analysis/Observed\\\_variable\\\_network\\\_recovery.png`, `Observed\\\_variable\\\_network\\\_analysis/Observed\\\_variable\\\_network\\\_metrics.xlsx`, `H1\\\_observed\\\_edge\\\_analysis/H1\\\_observed\\\_edge\\\_recovery\\\_and\\\_signal\\\_specificity\\\_1x3.png`, and `H1\\\_observed\\\_edge\\\_analysis/H1\\\_observed\\\_edge\\\_analysis.xlsx`.

The optional `Data\\\_generation.m` script recreates the simulated data design with 500 observed variables, five ground-truth hidden regulatory variables, one hidden global variable, and three repeats. The archived data should be used for exact reproduction of the reported results.

### 4.2 Paper astrocyte data experiment

Set the MATLAB current folder to `Paper\\\_astrocyte\\\_data\\\_method\\\_comparison`.

`Code\\\_Ca\\\_method\\\_comparison.m` requires the following files in that directory:

```text
S3\\\_Precision\\\_matrix.txt
S1\\\_whole\\\_gene\\\_list.xlsx
S5\\\_Ca\\\_associated\\\_reference\\\_label.xlsx
LVGGM\\\_L\\\_lowrank.txt
Data\\\_whole.txt
```

The script ranks the 2,958 non-prior genes using the HR-GGM observed--hidden-regulatory precision entries, a prior-selected LV-GGM component, and mean absolute correlation with the 42 prior genes. It produces top-K positive rates, enrichment folds, precision--recall results, permutation tests, bootstrap comparisons, and spectral diagnostics.

`Code\\\_statistic\\\_test.m` reproduces the threshold-based contingency analysis at an absolute HR-GGM score of 0.038 and reports the odds ratio, 95% confidence interval, Fisher's exact test, and chi-square test in the MATLAB command window.

The principal figure is `Method\\\_Comparison\\\_Output/MethodComparison\\\_Main\\\_ABC\\\_1x3.png`; ranking, resampling, enrichment, and spectral-diagnostic tables are written to the same output directory.

### 4.3 Supplementary analysis S1

**Robustness of observed-variable network and hidden-regulator association recovery.** Set the MATLAB current folder to `Supplementary\\\_simulation\\\_analysis\\\_S1` and run `Code\\\_simulation\\\_sensitivity.m`.

The experiment uses 500 observed variables and noise variance fractions of 0%, 10%, 25%, and 50%, with three repeats per condition. In the setting with hidden effects, the ground-truth network contains one hidden master regulatory variable, one hidden global variable, 600 undirected observed-variable edges, and 75 observed--hidden-regulatory edges. In the model-over-specification setting, data are generated without hidden effects, although HR-GGM is still fitted with both hidden variables. Observed-network recovery is therefore evaluated over 24 datasets, whereas observed--hidden-regulatory recovery is evaluated over the 12 datasets in which those effects are present.

The six-panel output reports observed-network F1 curves for both settings, threshold intervals retaining at least 90% of the maximum F1 score, observed--hidden-regulatory F1 curves, relative Frobenius error of the observed block, and sign-invariant relative Frobenius error of the observed--hidden-regulatory block. The principal output is `Sensitive\\\_O\\\_OH\\\_Analysis/Figure\\\_O\\\_OH\\\_combined\\\_reordered\\\_CleanLines.png`, corresponding to Supplementary Figure S1; numerical results are written to `Sensitive\\\_O\\\_OH\\\_Analysis/Threshold\\\_sensitive\\\_results.xlsx`.

### 4.4 Supplementary analysis S2

**Recovery of observed-variable edges under different sparsity levels.** Set the MATLAB current folder to `Supplementary\\\_simulation\\\_analysis\\\_S2` and run `Code\\\_HRGGM\\\_O\\\_sparsity.m`.

Each dataset contains 200 observed variables, one hidden master regulatory variable connected to 50 observed variables, and one hidden global variable connected to all 200 observed variables. The observed network retains 200 within-group edges and adds 0, 50, 200, or 1,000 between-group edges, giving 200, 250, 400, or 1,200 true undirected edges. These four network settings are crossed with noise variance fractions of 0%, 10%, and 25% and three repeats, giving 36 datasets of 10,000 samples each.

All possible observed-variable pairs are evaluated together. The three panels report normalized AUPRC, top-K recovery with K equal to the true number of undirected edges, and relative Frobenius error of the off-diagonal observed block. The principal output is `O\\\_Sparsity\\\_Performance\\\_Output/Figure\\\_O\\\_sparsity\\\_performance\\\_1x3.png`, corresponding to Supplementary Figure S2; condition-level and repeat-summary tables are written to the same output directory.

### 4.5 Supplementary analysis S3

**Validation of the data-driven parameter-selection procedure.** Set the MATLAB current folder to `Supplementary\\\_simulation\\\_analysis\\\_S3` and run `Code\\\_parameter\\\_selection\\\_validation.m`.

Each dataset contains 200 observed variables, one hidden master regulatory variable, one hidden global variable, and 10,000 samples. The ground truth contains 250 undirected observed-variable edges and 40 observed--hidden-regulatory edges; 16 of the 40 connected observed variables are supplied as priors and the remaining 24 are held out. Noise variance fractions of 0%, 20%, and 50% are evaluated over three matched repeats, with γ fixed at 0.9.

The standard HR-GGM procedure first selects λ by maximizing its stability criterion and α by minimizing its cross-validation criterion. Full-data λ and α scans are then evaluated only for post hoc comparison with ground-truth-based reference values; ground truth does not guide parameter selection or model fitting. The three panels show off-diagonal observed-block error across λ, the ground-truth-based score gap across α, and population Gaussian loss at the ground-truth score optimum versus the HR-GGM-selected α. The principal output is `Parameter\\\_Selection\\\_Validation\\\_OUTPUT/Figure\\\_Main\\\_ParameterSelection\\\_1x3.png`, corresponding to Supplementary Figure S3; scan and summary tables are written to the same output directory.

## 5\. Recommended execution order

For direct reproduction from the archived inputs and fitted matrices:

1. Merge `Data\\\_HRGGM\\\_paper` into `HRGGM\\\_paper` as described in Section 1.
2. Start MATLAB and set the current folder to one analysis module at a time.
3. Run the two main-text simulation analysis scripts.
4. Run `Code\\\_Ca\\\_method\\\_comparison.m`, followed by `Code\\\_statistic\\\_test.m`.
5. Run `Code\\\_simulation\\\_sensitivity.m`, `Code\\\_HRGGM\\\_O\\\_sparsity.m`, and `Code\\\_parameter\\\_selection\\\_validation.m` in the S1, S2, and S3 directories, respectively.

These modules do not need to be run in a single MATLAB session. Each primary analysis script writes to its own output directory.

To reproduce the complete computational pipeline from generated observations rather than from archived fitted matrices, run the relevant `Data\\\_generation.m`, fit HR-GGM with the module-specific `C\\\_code\\\_HRGGM`, run the required comparison methods, and then run the MATLAB analysis script. Read `C\\\_code\\\_HRGGM/README\\\_Data\\\_and\\\_Main\\\_mapping.txt` for the HR-GGM input-file format. Each module's `Main.cpp` is the authoritative run-specific configuration; do not replace it with a copy from another module.

## 6\. Important reproducibility notes

* Keep all matrix row and column orders unchanged.
* C++ index files use zero-based indices; the MATLAB analysis scripts convert them to one-based indices where required.
* `Theta\\\_full.txt` is simulation ground truth and is not used by HR-GGM parameter selection or fitting.
* The archived fitted matrices and intermediate score files are the inputs for exact figure reproduction.
* Several scripts overwrite files in their output directories. `Code\\\_parameter\\\_selection\\\_validation.m` recreates its output directory, and the simulation generators have `reset\\\_output\\\_folder = true`. Preserve a copy of archived outputs before intentionally regenerating data or results.
* Do not rename source dataset directories or legacy table columns solely to match manuscript terminology; the analysis scripts already translate exported descriptions where appropriate.

## 7\. Citation and third-party software

If you use this code or the accompanying data, please cite the HR-GGM paper:

> \\\*\\\*HR-GGM paper citation to be added after publication.\\\*\\\*

The comparison directory contains code or wrappers associated with LV-GGM, TGL, GRNBoost2, SCING, and RegDiffusion. Please also cite the corresponding original methods when using or rerunning them. The bundled Eigen and LogdetPPA components retain their own license and citation files where provided.

## 8\. Data and code availability

* Code: this GitHub repository (`HRGGM\\\_paper`)
* Data and archived fitted results: [`Data\\\_HRGGM\\\_paper.zip` on Zenodo](https://zenodo.org/records/22217595?preview=1&token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjMwZjVkZTQ3LTAxMzYtNDEyNy04ODllLTUzZTRkNjQ5ZTU2MiIsImRhdGEiOnt9LCJyYW5kb20iOiI5NmU1MTNhMGQyNjA3ZTNlY2NhNzhkNTgzYzliNmQ0MSJ9.4h8NHWIfglk8OPAE7i7CWuAb8tDwtuTK2Bnncz2Y7SdJj_6rkwYNO_lcI6AqRKLJzzje_dmQGqy6YtR--EteSg) (draft preview)

