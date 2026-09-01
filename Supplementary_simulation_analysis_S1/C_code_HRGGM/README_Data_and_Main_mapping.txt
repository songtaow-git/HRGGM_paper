HR-GGM C++ Project
Data Directory Layout and Main.cpp Input/Output Mapping
=======================================================

Terminology Note
----------------

The C++ implementation retains legacy identifiers, and these identifiers
should not be renamed. Their correspondence to the current manuscript is:

- observed variables (`O`) correspond to `local`, `L`, and `Nlocal` in code;
- hidden master regulatory variables (`H`) correspond to `hidden`, `H`, and
  `Nhidden` in code;
- hidden global variables (`G`) correspond to `global`, `G`, and `Nglobal` in
  code.

Literal C++ variable names, function names, directory names, and output file
names are reproduced exactly below. The surrounding explanations use the
current manuscript terminology.

1. Project Layout
-----------------

After building the project, place the complete `Data` directory inside the
compiled project directory.

Recommended layout:

Project/
└── build/
    ├── C_GRN                  Linux executable
    ├── C_GRN.exe              Windows executable name may differ
    ├── Data/
    │   ├── Alpha_list.txt
    │   ├── Data_whole.bin
    │   ├── Idx_g_list.txt
    │   ├── Idx_h_list.txt
    │   ├── Lambda_list.txt
    │   └── Theta_full.txt        Used only for simulation evaluation
    └── Result/                Created automatically at runtime

If the executable is located in a subdirectory, for example:

Project/build/Release/C_GRN.exe

the `Data` directory may remain at:

Project/build/Data/

`Main.cpp` starts from the executable directory and searches upward for the
nearest directory containing `Data`. Once found, the program creates `Result`
beside `Data`.

For example:

Project/build/
├── Data/
├── Result/
└── Release/
    └── C_GRN.exe


2. Path Resolution in Main.cpp
------------------------------

The executable location is determined by:

    std::filesystem::path exe_path = std::filesystem::absolute(argv[0]);
    std::filesystem::path exe_dir = exe_path.parent_path();

The code then searches upward until it finds:

    base_dir / "Data"

The main directories are defined as:

    data_dir = base_dir / "Data";
    out_dir  = base_dir / "Result";

Two subdirectories are created automatically:

    Result/lambda_ratio/
    Result/Alpha_metric/


3. Data Files and Their Roles in Main.cpp
-----------------------------------------

3.1 Data_whole.bin
~~~~~~~~~~~~~~~~~~

Main.cpp variables:

    Path_data
    Data
    Nlocal
    Nfeature

Read by:

    std::string Path_data = (data_dir / "Data_whole.bin").string();
    RowMat Data = Load_bin_f64_data(Path_data);

Purpose:

- Stores the observed-variable data matrix.
- Matrix orientation is observed variables × samples.
- `Data.rows()` determines `Nlocal`, the number of observed variables.
- The data are used to construct cross-validation folds.
- The data are used to compute fold-specific extended data matrices
  (`LHG_train` and `LHG_test` in code) and covariance matrices.
- The full dataset is used for the final model fit.

Main functions using this input:

    MakeKfoldplan_reorderdata
    Precompute_foldzscores
    Compute_fold_lhg_cov_parallel
    Solution_theta_and_hg_factor


3.2 Lambda_list.txt
~~~~~~~~~~~~~~~~~~~

Main.cpp variables:

    Path_lambda
    Lambda_list
    Lambda_fix
    Lambda_old

Read by:

    std::string Path_lambda = (data_dir / "Lambda_list.txt").string();
    std::vector<double> Lambda_list = ReadParameterFromTxt(Path_lambda);

Purpose:

- Stores candidate lambda values.
- Used in the initial lambda search.
- Used again during alternating alpha/lambda optimization.
- The selected value is stored in `Lambda_fix`.

Main functions using this input:

    Compute_orig_lambda_parallel
    Score_orig.Compute_score
    Compute_fold_lambda_parallel
    Lam_metric.Compute_score

Related output directory:

    Result/lambda_ratio/


3.3 Alpha_list.txt
~~~~~~~~~~~~~~~~~~

Main.cpp variables:

    Path_alpha
    Alpha_list
    Alpha_fix
    Alpha_old
    Alpha_retention_tau

Read by:

    std::string Path_alpha = (data_dir / "Alpha_list.txt").string();
    std::vector<double> Alpha_list = ReadParameterFromTxt(Path_alpha);

Purpose:

- Stores candidate alpha values.
- Evaluated while the current lambda value is fixed.
- `Alpha_retention_tau` is used in the alpha-selection score.
- The selected value is stored in `Alpha_fix`.

Main functions using this input:

    Compute_fold_alpha_parallel
    Alp_metric.Compute_score

Related output directory:

    Result/Alpha_metric/


3.4 Idx_h_list.txt
~~~~~~~~~~~~~~~~~~

Main.cpp variables:

    Path_idx_h
    valid_idx_h

Read by:

    std::string Path_idx_h = (data_dir / "Idx_h_list.txt").string();
    std::vector<int> valid_idx_h = ReadidxFromTxt(Path_idx_h);

Purpose:

- Stores indices of prior observed variables used to initialize the hidden
  master regulatory variable or variables.
- Indices must use zero-based C++ indexing.
- Used to initialize the hidden master regulatory variable or variables within
  the joint hidden/global construction procedure.
- Used during alpha selection, lambda selection, and final model fitting.

Main functions using this input:

    Compute_fold_lhg_cov_parallel
    Compute_fold_alpha_parallel
    Compute_fold_lambda_parallel
    Solution_theta_and_hg_factor


3.5 Idx_g_list.txt
~~~~~~~~~~~~~~~~~~

Main.cpp variables:

    Path_idx_g
    valid_idx_g

Read by:

    std::string Path_idx_g = (data_dir / "Idx_g_list.txt").string();
    std::vector<int> valid_idx_g = ReadidxFromTxt(Path_idx_g);

Purpose:

- Stores indices of observed variables used to initialize the hidden global
  variable or variables.
- Indices must use zero-based C++ indexing.
- Used to initialize the hidden global variable or variables within the joint
  hidden/global construction procedure.
- Used during alpha selection, lambda selection, and final model fitting.

Main functions using this input:

    Compute_fold_lhg_cov_parallel
    Compute_fold_alpha_parallel
    Compute_fold_lambda_parallel
    Solution_theta_and_hg_factor


3.6 Theta_full.txt
~~~~~~~~~~~~~~~~~~

`Theta_full.txt` is not read by `Main.cpp`.

Purpose:

- Stores the ground-truth precision matrix for simulated data.
- Used only for downstream evaluation of `Result/Theta_final.txt`.
- Does not affect parameter selection or model fitting.

Typical dimensions:

- Simulation containing observed, hidden regulatory, and hidden global
  variables:
      (Nlocal + Nhidden + Nglobal) ×
      (Nlocal + Nhidden + Nglobal)

- Simulation containing observed variables only:
      Nlocal × Nlocal


4. User-Adjustable Parameters in Main.cpp
-----------------------------------------

    Max_iter
        Maximum number of iterations used by each precision-matrix
        optimization step.

    Converge_thre
        Convergence tolerance for numerical optimization.

    N_thread
        Number of parallel threads used during parameter selection.

    Blockcols
        Number of sample columns processed per covariance-computation block.

    Fold_k
        Number of cross-validation folds.

    Nhidden
        Number of hidden master regulatory variables included in the model.

    Nglobal
        Number of hidden global variables included in the model.

    Gamma
        Mixing weight for regularization of the observed--hidden regulatory
        block. The C++ variable name remains `Gamma`.

    Alpha_retention_tau
        Retention threshold used in the alpha-selection score.


5. Output Files
---------------

The program creates:

Result/
├── Theta_final.txt
├── lambda_ratio/
└── Alpha_metric/


5.1 Theta_final.txt
~~~~~~~~~~~~~~~~~~~

Generated by:

    Solution_theta_and_hg_factor(...)
    WriteSparseForMatlab(Input_whole.Get_theta(), ...)

Contents:

- Final estimated extended precision matrix.
- Its dimension is `Nfeature × Nfeature`, where
  `Nfeature = Nlocal + Nhidden + Nglobal`.
- Rows and columns are ordered as observed variables, hidden master regulatory
  variables, and hidden global variables.
- Computed from the complete input dataset.
- Uses the selected `Lambda_fix`, `Alpha_fix`, and `Gamma` values.


5.2 Result/lambda_ratio/
~~~~~~~~~~~~~~~~~~~~~~~~

Typical files:

    Score_norm_0.txt
    Score_ratio_0.txt
    Score_norm<loop>.txt
    Score_ratio<loop>.txt

Interpretation:

- Files ending in `_0` correspond to the initial lambda search.
- Files indexed by `loop` correspond to lambda evaluation during alternating
  alpha/lambda optimization.


5.3 Result/Alpha_metric/
~~~~~~~~~~~~~~~~~~~~~~~~

Typical files generated during each alpha-search iteration:

    Likelihood_mat_<loop>.txt
    Entropy_mat_<loop>.txt
    Global_magnitude_mat_<loop>.txt
    Retention_vec_<loop>.txt
    Adjusted_entropy_vec_<loop>.txt
    Score_vec_<loop>.txt

These files store the corresponding matrices and vectors from the
`Alpha_metric` object for each alpha-selection iteration.


6. Pre-Run Checklist
--------------------

Confirm that the `Data` directory contains exactly:

    Alpha_list.txt
    Data_whole.bin
    Idx_g_list.txt
    Idx_h_list.txt
    Lambda_list.txt
    Theta_full.txt

`Main.cpp` directly reads the first five files.

`Theta_full.txt` is retained for downstream ground-truth evaluation only.

Do not rename these files or change their letter case. On Windows, file
extensions may be hidden in File Explorer even though the files are stored
with the `.txt` extension.
