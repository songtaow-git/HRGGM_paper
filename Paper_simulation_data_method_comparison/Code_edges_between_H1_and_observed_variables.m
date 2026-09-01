%% Recovery of edges between H1 and observed variables
% Evaluates observed-variable network recovery, held-out H1 target recovery, and
% H1-related signal composition across three simulation repeats.
%
% Network dimensions:
%   500 observed variables
%   5 ground-truth hidden regulatory variables, with H1 as the target
%   1 ground-truth hidden global variable
%
% Matrix dimensions:
%   Ground truth: 506 x 506
%   HR-GGM:       502 x 502
%   LV-GGM:       500 x 500 observed-variable precision matrix plus latent input
%
% The H1 design contains 75 associated observed variables, including 30 prior
% variables and 45 held-out targets. Prior indices are read as C++ 0-based
% indices and converted to MATLAB 1-based indices.
%
% Outputs:
%   H1_observed_edge_recovery_and_signal_specificity_1x3.png/.pdf
%   H1_observed_edge_analysis.xlsx

clear;
clc;
close all;

%% ============================== Settings ================================

Data_folder  = fullfile('Data');
GT_folder    = fullfile('GT');
HR_folder    = fullfile('HRGGM');
LV_folder    = fullfile('LVGGM');
Prior_folder = fullfile('Prior_genes');

output_dir = fullfile('H1_observed_edge_analysis');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

Data_prefix       = 'Data_whole';
GT_prefix         = 'GT';
HR_prefix         = 'HRGGM';
LV_sparse_prefix  = 'LVGGM';
LV_latent_prefix  = 'LVGGM_latent';
Prior_prefix      = 'Idx_h_list';

Nlocal    = 500;
NhiddenGT = 5;
NglobalGT = 1;
NhiddenHR = 1;
NglobalHR = 1;
Nrep      = 3;

GT_expected_dim = Nlocal + NhiddenGT + NglobalGT;   % 506
HR_expected_dim = Nlocal + NhiddenHR + NglobalHR;   % 502
LV_expected_dim = Nlocal;                           % 500

% GT indexing
idxL_GT  = 1:Nlocal;
idxH_GT  = (Nlocal + 1):(Nlocal + NhiddenGT);        % 501:505
idxH1_GT = idxH_GT(1);                              % 501
idxG_GT  = GT_expected_dim;                         % 506

% HR-GGM indexing
idxL_HR  = 1:Nlocal;
idxH1_HR = Nlocal + 1;                              % 501

% Prior index count.
expected_prior_count        = 30;

GT_zero_tolerance = 1e-12;
target_recall     = 0.95;

candidate_multipliers_LL = [0.5, 1, 2, 3, 5, 10, 20];
candidate_multipliers_H1 = [0.5, 1, 2, 3, 5, 10];
recall_grid = (0:0.001:1)';

% Canonical LV-GGM latent-component handling
relative_eigenvalue_tolerance = 1e-8;
absolute_eigenvalue_tolerance = 1e-12;
square_symmetry_tolerance = 1e-5;
auto_orient_symmetric_latent_matrix = true;

% Figure settings
font_name          = 'Arial';
font_size_axis     = 8.5;
font_size_title    = 9.5;
font_size_legend   = 7.5;
font_size_annotation = 7.2;
figure_resolution  = 1200;
figure2_width_cm   = 17.8;   % 178 mm, Bioinformatics double-column width
figure2_height_cm  = 8.2;

% Methods evaluated for observed-variable network recovery.
LL_method_names = {'HR-GGM', 'LV-GGM'};
Nmethod_LL = numel(LL_method_names);

LL_method_colors = [ ...
    0.0000, 0.4470, 0.7410; ...
    0.8500, 0.3250, 0.0980];

LL_method_markers = {'o', 's'};

% Methods evaluated for held-out H1 recovery.
H1_method_names = { ...
    'HR-GGM', ...
    'LV-GGM', ...
    'Prior-correlation'};
Nmethod_H1 = numel(H1_method_names);

H1_method_colors = [ ...
    0.0000, 0.4470, 0.7410; ...
    0.8500, 0.3250, 0.0980; ...
    0.4660, 0.6740, 0.1880];

H1_method_markers = {'o', 's', '^'};

% Estimated H1-related signals used in the composition analysis.
% The composition analysis includes the HR-GGM H1 and LV-GGM H1 estimates.
estimate_names = { ...
    'HR-GGM hidden signal'; ...
    'LV-GGM oracle H1 component'};

% Ground-truth reference signals used in the composition analysis.
reference_names = { ...
    'GT target H1 signal'; ...
    'GT other-hidden signal'; ...
    'GT global signal'};

Nestimate  = numel(estimate_names);
Nreference = numel(reference_names);

%% ========================= Allocate results =============================

% Observed-variable network truth information by repeat.
LL_TrueEdges          = nan(Nrep, 1);
LL_TotalPossibleEdges = nan(Nrep, 1);
LL_PositiveFraction   = nan(Nrep, 1);

% Observed-variable network metrics: method by repeat.
LL_AUPRC            = nan(Nmethod_LL, Nrep);
LL_NormalizedAUPRC  = nan(Nmethod_LL, Nrep);
LL_TopKRecovery     = nan(Nmethod_LL, Nrep);
LL_Budget95Edges    = nan(Nmethod_LL, Nrep);
LL_Budget95Factor   = nan(Nmethod_LL, Nrep);
LL_PrecisionAt95    = nan(Nmethod_LL, Nrep);
LL_FalsePositive95  = nan(Nmethod_LL, Nrep);
LL_PR_grid = nan(numel(recall_grid), Nmethod_LL, Nrep);
LL_candidate_recovery = nan( ...
    numel(candidate_multipliers_LL), Nmethod_LL, Nrep);

% Held-out H1 truth information by repeat
H1_TotalTargets           = nan(Nrep, 1);
H1_PriorCount             = nan(Nrep, 1);
H1_HeldoutTrueEdges       = nan(Nrep, 1);
H1_TotalEvaluatedFeatures = nan(Nrep, 1);
H1_PositiveFraction       = nan(Nrep, 1);

% Held-out H1 metrics: H1 method x repeat
H1_AUPRC              = nan(Nmethod_H1, Nrep);
H1_NormalizedAUPRC    = nan(Nmethod_H1, Nrep);
H1_TopKRecovery       = nan(Nmethod_H1, Nrep);
H1_Budget95Edges      = nan(Nmethod_H1, Nrep);
H1_Budget95Factor     = nan(Nmethod_H1, Nrep);
H1_PrecisionAt95      = nan(Nmethod_H1, Nrep);
H1_FalsePositive95    = nan(Nmethod_H1, Nrep);
H1_PearsonWeightCorr  = nan(Nmethod_H1, Nrep);
H1_SpearmanWeightCorr = nan(Nmethod_H1, Nrep);
H1_PR_grid = nan(numel(recall_grid), Nmethod_H1, Nrep);
H1_candidate_recovery = nan( ...
    numel(candidate_multipliers_H1), Nmethod_H1, Nrep);

% Signal comparison: repeat x estimate x reference
SignalCorrelation = nan(Nrep, Nestimate, Nreference);
SignalBeta        = nan(Nrep, Nestimate, Nreference);
SignalProportion  = nan(Nrep, Nestimate, Nreference);
SignalR2          = nan(Nrep, Nestimate);

% LV-GGM latent selection diagnostics
LV_RetainedRank             = nan(Nrep, 1);
LV_OracleSelectedComponent  = nan(Nrep, 1);
LV_OracleCorrelationH1      = nan(Nrep, 1);
LV_PriorSelectedComponent   = nan(Nrep, 1);
LV_PriorEnrichment          = nan(Nrep, 1);

% GT structural diagnostics
GT_HiddenBlockOffdiagRatio  = nan(Nrep, 1);
GT_HiddenGlobalCouplingNorm = nan(Nrep, 1);

%% =========================== Main analysis ==============================

for rep = 1:Nrep

    %% -------------------------- Resolve files ---------------------------

    Data_file = resolve_data_file(Data_folder, Data_prefix, rep);
    GT_file = resolve_data_file(GT_folder, GT_prefix, rep);
    HR_file = resolve_data_file(HR_folder, HR_prefix, rep);
    LV_sparse_file = resolve_data_file( ...
        LV_folder, LV_sparse_prefix, rep);
    LV_latent_file = resolve_data_file( ...
        LV_folder, LV_latent_prefix, rep);
    Prior_file = resolve_data_file( ...
        Prior_folder, Prior_prefix, rep);

    %% --------------------------- Read files -----------------------------

    Data_raw = read_numeric_array(Data_file);
    Theta_GT = read_numeric_array(GT_file);
    Theta_HR = read_numeric_array(HR_file);
    Theta_LV = read_numeric_array(LV_sparse_file);
    LV_latent_raw = read_numeric_array(LV_latent_file);

    assert_matrix_size(Theta_GT, GT_expected_dim, GT_expected_dim, ...
        sprintf('GT repeat %d', rep));
    assert_matrix_size(Theta_HR, HR_expected_dim, HR_expected_dim, ...
        sprintf('HR-GGM repeat %d', rep));
    assert_matrix_size(Theta_LV, LV_expected_dim, LV_expected_dim, ...
        sprintf('LV-GGM sparse repeat %d', rep));

    Theta_GT = symmetrize_matrix(Theta_GT);
    Theta_HR = symmetrize_matrix(Theta_HR);
    Theta_LV = symmetrize_matrix(Theta_LV);

    X_local = prepare_local_feature_data( ...
        Data_raw, ...
        Nlocal, ...
        [Nlocal, HR_expected_dim, GT_expected_dim], ...
        rep);

    %% -------------------- Read and convert priors -----------------------

    prior_idx_cpp = read_cpp_index_vector(Prior_file);

    if numel(prior_idx_cpp) ~= expected_prior_count
        error([ ...
            'Repeat %d prior file must contain exactly %d indices, ', ...
            'but %d were detected.'], ...
            rep, expected_prior_count, numel(prior_idx_cpp));
    end

    if any(prior_idx_cpp < 0) || any(prior_idx_cpp > Nlocal - 1)
        error([ ...
            'Repeat %d contains invalid C++ prior indices. ', ...
            'Allowed range is 0 to %d.'], ...
            rep, Nlocal - 1);
    end

    if numel(unique(prior_idx_cpp)) ~= expected_prior_count
        error('Duplicate prior indices were found in repeat %d.', rep);
    end

    % C++ 0-based -> MATLAB 1-based
    prior_idx_matlab = sort(prior_idx_cpp(:) + 1);
    nonprior_idx_matlab = setdiff((1:Nlocal)', prior_idx_matlab);

    %% --------------------------------------------------------------------
    % Observed-variable network recovery
    %% --------------------------------------------------------------------

    GT_LL = max( ...
        abs(Theta_GT(idxL_GT, idxL_GT)), ...
        abs(Theta_GT(idxL_GT, idxL_GT)'));

    HR_LL = max( ...
        abs(Theta_HR(idxL_HR, idxL_HR)), ...
        abs(Theta_HR(idxL_HR, idxL_HR)'));

    LV_LL = max(abs(Theta_LV), abs(Theta_LV'));

    GT_LL(1:Nlocal+1:end) = 0;
    HR_LL(1:Nlocal+1:end) = 0;
    LV_LL(1:Nlocal+1:end) = 0;

    upper_mask = triu(true(Nlocal), 1);

    LL_truth = GT_LL(upper_mask) > GT_zero_tolerance;
    LL_truth = logical(LL_truth(:));

    K_LL = sum(LL_truth);
    Npairs_LL = numel(LL_truth);
    prevalence_LL = K_LL / Npairs_LL;

    if K_LL == 0
        error('No true LL edges were found in GT repeat %d.', rep);
    end

    LL_TrueEdges(rep)          = K_LL;
    LL_TotalPossibleEdges(rep) = Npairs_LL;
    LL_PositiveFraction(rep)   = prevalence_LL;

    LL_score_cells = { ...
        HR_LL(upper_mask); ...
        LV_LL(upper_mask)};

    for method_id = 1:Nmethod_LL

        current_scores = LL_score_cells{method_id};
        current_scores = current_scores(:);

        Ranking = calculate_ranking_metrics( ...
            current_scores, ...
            LL_truth, ...
            K_LL, ...
            target_recall, ...
            candidate_multipliers_LL, ...
            recall_grid);

        LL_AUPRC(method_id, rep) = Ranking.AUPRC;
        LL_NormalizedAUPRC(method_id, rep) = ...
            (Ranking.AUPRC - prevalence_LL) / ...
            max(1 - prevalence_LL, eps);
        LL_TopKRecovery(method_id, rep) = Ranking.TopKRecovery;
        LL_Budget95Edges(method_id, rep) = Ranking.TargetBudgetEdges;
        LL_Budget95Factor(method_id, rep) = Ranking.TargetBudgetFactor;
        LL_PrecisionAt95(method_id, rep) = Ranking.PrecisionAtTarget;
        LL_FalsePositive95(method_id, rep) = Ranking.FalsePositiveBurden;
        LL_PR_grid(:, method_id, rep) = Ranking.PrecisionOnRecallGrid;
        LL_candidate_recovery(:, method_id, rep) = ...
            Ranking.CandidateRecovery;

    end

    %% --------------------------------------------------------------------
    % Canonical LV-GGM latent components
    %% --------------------------------------------------------------------

    LV_components = prepare_canonical_components( ...
        LV_latent_raw, ...
        Nlocal, ...
        relative_eigenvalue_tolerance, ...
        absolute_eigenvalue_tolerance, ...
        square_symmetry_tolerance, ...
        auto_orient_symmetric_latent_matrix, ...
        'LV-GGM', ...
        rep);

    LV_RetainedRank(rep) = LV_components.n_components;

    %% --------------------------------------------------------------------
    % Ground-truth H1 target categories
    %% --------------------------------------------------------------------

    GT_H1_edge_scores = abs(Theta_GT(idxL_GT, idxH1_GT));
    GT_H1_edge_scores = GT_H1_edge_scores(:);

    H1_target_mask = GT_H1_edge_scores > GT_zero_tolerance;
    H1_target_mask = logical(H1_target_mask(:));

    if any(~H1_target_mask(prior_idx_matlab))
        bad_prior = prior_idx_matlab(~H1_target_mask(prior_idx_matlab));
        error([ ...
            'Repeat %d contains prior indices that are not true H1 ', ...
            'targets. MATLAB indices: %s'], ...
            rep, mat2str(bad_prior(:)'));
    end

    total_H1_targets = sum(H1_target_mask);

    heldout_H1_mask = H1_target_mask;
    heldout_H1_mask(prior_idx_matlab) = false;

    eval_truth = heldout_H1_mask(nonprior_idx_matlab);
    eval_truth = logical(eval_truth(:));

    K_H1 = sum(eval_truth);
    N_eval_H1 = numel(eval_truth);
    prevalence_H1 = K_H1 / N_eval_H1;
    if K_H1 == 0
        error('No held-out H1 targets were found in repeat %d.', rep);
    end

    H1_TotalTargets(rep) = total_H1_targets;
    H1_PriorCount(rep) = numel(prior_idx_matlab);
    H1_HeldoutTrueEdges(rep) = K_H1;
    H1_TotalEvaluatedFeatures(rep) = N_eval_H1;
    H1_PositiveFraction(rep) = prevalence_H1;

    %% --------------------------------------------------------------------
    % Prior-matched recovery of held-out H1-associated variables
    %% --------------------------------------------------------------------

    % HR-GGM scores from observed-H1 precision entries.
    HR_H1_scores_full = abs(Theta_HR(idxL_HR, idxH1_HR));
    HR_H1_scores_full = HR_H1_scores_full(:);

    % Select the LV-GGM canonical component using the prior set.
    LV_prior = select_component_by_prior( ...
        LV_components, ...
        prior_idx_matlab);

    LV_PriorSelectedComponent(rep) = LV_prior.selected_rank;
    LV_PriorEnrichment(rep) = LV_prior.selected_prior_enrichment;

    LV_H1_scores_full = LV_prior.selected_feature_score(:);
    % Mean absolute Pearson correlation with the 30 prior variables.
    PriorCorr_H1_scores_full = calculate_prior_correlation_scores( ...
        X_local, ...
        prior_idx_matlab);

    H1_score_cells = { ...
        HR_H1_scores_full; ...
        LV_H1_scores_full; ...
        PriorCorr_H1_scores_full};

    GT_H1_eval_weights = GT_H1_edge_scores(nonprior_idx_matlab);

    for method_id = 1:Nmethod_H1

        current_full_scores = H1_score_cells{method_id};
        current_eval_scores = current_full_scores(nonprior_idx_matlab);

        Ranking = calculate_ranking_metrics( ...
            current_eval_scores, ...
            eval_truth, ...
            K_H1, ...
            target_recall, ...
            candidate_multipliers_H1, ...
            recall_grid);

        H1_AUPRC(method_id, rep) = Ranking.AUPRC;
        H1_NormalizedAUPRC(method_id, rep) = ...
            (Ranking.AUPRC - prevalence_H1) / ...
            max(1 - prevalence_H1, eps);
        H1_TopKRecovery(method_id, rep) = Ranking.TopKRecovery;
        H1_Budget95Edges(method_id, rep) = Ranking.TargetBudgetEdges;
        H1_Budget95Factor(method_id, rep) = Ranking.TargetBudgetFactor;
        H1_PrecisionAt95(method_id, rep) = Ranking.PrecisionAtTarget;
        H1_FalsePositive95(method_id, rep) = Ranking.FalsePositiveBurden;
        H1_PR_grid(:, method_id, rep) = Ranking.PrecisionOnRecallGrid;
        H1_candidate_recovery(:, method_id, rep) = ...
            Ranking.CandidateRecovery;

        H1_PearsonWeightCorr(method_id, rep) = safe_corr( ...
            current_eval_scores, GT_H1_eval_weights, 'Pearson');

        H1_SpearmanWeightCorr(method_id, rep) = safe_corr( ...
            current_eval_scores, GT_H1_eval_weights, 'Spearman');

    end

    %% --------------------------------------------------------------------
    % Hidden regulatory and global observed-space signal specificity
    %
    % Construct each ground-truth hidden signal separately.
    % Each signal uses the corresponding precision-matrix blocks.
    %% --------------------------------------------------------------------

    GT_hidden_components = cell(NhiddenGT, 1);

    for h = 1:NhiddenGT
        hidden_index = idxH_GT(h);
        GT_hidden_components{h} = latent_component_signal( ...
            Theta_GT, idxL_GT, hidden_index);
    end

    S_GT_H1 = GT_hidden_components{1};

    S_GT_other_hidden = zeros(Nlocal, Nlocal);
    for h = 2:NhiddenGT
        S_GT_other_hidden = ...
            S_GT_other_hidden + GT_hidden_components{h};
    end
    S_GT_other_hidden = symmetrize_matrix(S_GT_other_hidden);

    S_GT_global = latent_component_signal(Theta_GT, idxL_GT, idxG_GT);
    S_HR_H1     = latent_component_signal(Theta_HR, idxL_HR, idxH1_HR);

    % Select the LV-GGM component with maximum correlation to the GT H1 signal.
    LV_oracle = select_oracle_component(LV_components, S_GT_H1);

    LV_OracleSelectedComponent(rep) = LV_oracle.selected_rank;
    LV_OracleCorrelationH1(rep) = LV_oracle.selected_correlation_H1;

    S_LV_oracle_H1 = LV_oracle.selected_component;

    GT_HH = Theta_GT(idxH_GT, idxH_GT);
    GT_HH_offdiag = GT_HH - diag(diag(GT_HH));
    GT_HiddenBlockOffdiagRatio(rep) = ...
        norm(GT_HH_offdiag, 'fro') / max(norm(GT_HH, 'fro'), eps);

    GT_HiddenGlobalCouplingNorm(rep) = ...
        norm(Theta_GT(idxH_GT, idxG_GT), 'fro');

    estimated_signals = { ...
        S_HR_H1; ...
        S_LV_oracle_H1};

    reference_signals = { ...
        S_GT_H1; ...
        S_GT_other_hidden; ...
        S_GT_global};

    for e = 1:Nestimate
        for r = 1:Nreference
            SignalCorrelation(rep, e, r) = ...
                upper_triangle_correlation( ...
                    estimated_signals{e}, ...
                    reference_signals{r});
        end
    end

    estimated_row_scores = cell(Nestimate, 1);
    for e = 1:Nestimate
        estimated_row_scores{e} = ...
            sum(abs(estimated_signals{e}), 2);
    end

    reference_row_scores = [ ...
        sum(abs(S_GT_H1), 2), ...
        sum(abs(S_GT_other_hidden), 2), ...
        sum(abs(S_GT_global), 2)];

    for e = 1:Nestimate
        regression_result = standardized_component_regression( ...
            estimated_row_scores{e}, ...
            reference_row_scores);

        SignalBeta(rep, e, :) = regression_result.beta;
        SignalProportion(rep, e, :) = regression_result.proportion;
        SignalR2(rep, e) = regression_result.R2;
    end
end

%% ========================== Output tables ===============================

% ---------------- Observed-variable network results by repeat ------------

Nrow_LL = Nmethod_LL * Nrep;
LL_Table_Repeat = zeros(Nrow_LL, 1);
LL_Table_Method = strings(Nrow_LL, 1);
LL_Table_TrueEdges = zeros(Nrow_LL, 1);
LL_Table_TotalPairs = zeros(Nrow_LL, 1);
LL_Table_Prevalence = zeros(Nrow_LL, 1);
LL_Table_AUPRC = zeros(Nrow_LL, 1);
LL_Table_NormalizedAUPRC = zeros(Nrow_LL, 1);
LL_Table_TopK = zeros(Nrow_LL, 1);
LL_Table_Budget95Edges = zeros(Nrow_LL, 1);
LL_Table_Budget95Factor = zeros(Nrow_LL, 1);
LL_Table_Precision95 = zeros(Nrow_LL, 1);
LL_Table_FP95 = zeros(Nrow_LL, 1);

row_id = 0;
for rep = 1:Nrep
    for method_id = 1:Nmethod_LL
        row_id = row_id + 1;
        LL_Table_Repeat(row_id) = rep;
        LL_Table_Method(row_id) = LL_method_names{method_id};
        LL_Table_TrueEdges(row_id) = LL_TrueEdges(rep);
        LL_Table_TotalPairs(row_id) = LL_TotalPossibleEdges(rep);
        LL_Table_Prevalence(row_id) = LL_PositiveFraction(rep);
        LL_Table_AUPRC(row_id) = LL_AUPRC(method_id, rep);
        LL_Table_NormalizedAUPRC(row_id) = ...
            LL_NormalizedAUPRC(method_id, rep);
        LL_Table_TopK(row_id) = LL_TopKRecovery(method_id, rep);
        LL_Table_Budget95Edges(row_id) = LL_Budget95Edges(method_id, rep);
        LL_Table_Budget95Factor(row_id) = LL_Budget95Factor(method_id, rep);
        LL_Table_Precision95(row_id) = LL_PrecisionAt95(method_id, rep);
        LL_Table_FP95(row_id) = LL_FalsePositive95(method_id, rep);
    end
end

LL_PerRepeat = table( ...
    LL_Table_Repeat, LL_Table_Method, ...
    LL_Table_TrueEdges, LL_Table_TotalPairs, LL_Table_Prevalence, ...
    LL_Table_AUPRC, LL_Table_NormalizedAUPRC, LL_Table_TopK, ...
    LL_Table_Budget95Edges, LL_Table_Budget95Factor, ...
    LL_Table_Precision95, LL_Table_FP95, ...
    'VariableNames', { ...
        'Repeat', 'Method', 'TrueEdges', 'TotalPossibleEdges', ...
        'PositiveFraction', 'AUPRC', 'NormalizedAUPRC', ...
        'TopKRecovery', 'Budget95Edges', 'Budget95Factor', ...
        'PrecisionAt95Recall', 'FalsePositiveBurden95'});

LL_Summary = make_method_summary_table( ...
    LL_method_names, ...
    LL_AUPRC, ...
    LL_NormalizedAUPRC, ...
    LL_TopKRecovery, ...
    LL_Budget95Factor, ...
    LL_PrecisionAt95, ...
    LL_FalsePositive95, ...
    [], []);

LL_CandidateRecovery = make_candidate_table( ...
    LL_method_names, candidate_multipliers_LL, LL_candidate_recovery);

LL_MeanPRCurve = make_pr_curve_table( ...
    LL_method_names, recall_grid, LL_PR_grid);

% ------------------------- H1 per repeat ---------------------------------

Nrow_H1 = Nmethod_H1 * Nrep;
H1_Table_Repeat = zeros(Nrow_H1, 1);
H1_Table_Method = strings(Nrow_H1, 1);
H1_Table_TotalTargets = zeros(Nrow_H1, 1);
H1_Table_Priors = zeros(Nrow_H1, 1);
H1_Table_HeldoutTargets = zeros(Nrow_H1, 1);
H1_Table_Evaluated = zeros(Nrow_H1, 1);
H1_Table_Prevalence = zeros(Nrow_H1, 1);
H1_Table_AUPRC = zeros(Nrow_H1, 1);
H1_Table_NormalizedAUPRC = zeros(Nrow_H1, 1);
H1_Table_TopK = zeros(Nrow_H1, 1);
H1_Table_Budget95Edges = zeros(Nrow_H1, 1);
H1_Table_Budget95Factor = zeros(Nrow_H1, 1);
H1_Table_Precision95 = zeros(Nrow_H1, 1);
H1_Table_FP95 = zeros(Nrow_H1, 1);
H1_Table_Pearson = zeros(Nrow_H1, 1);
H1_Table_Spearman = zeros(Nrow_H1, 1);

row_id = 0;
for rep = 1:Nrep
    for method_id = 1:Nmethod_H1
        row_id = row_id + 1;
        H1_Table_Repeat(row_id) = rep;
        H1_Table_Method(row_id) = H1_method_names{method_id};
        H1_Table_TotalTargets(row_id) = H1_TotalTargets(rep);
        H1_Table_Priors(row_id) = H1_PriorCount(rep);
        H1_Table_HeldoutTargets(row_id) = H1_HeldoutTrueEdges(rep);
        H1_Table_Evaluated(row_id) = H1_TotalEvaluatedFeatures(rep);
        H1_Table_Prevalence(row_id) = H1_PositiveFraction(rep);
        H1_Table_AUPRC(row_id) = H1_AUPRC(method_id, rep);
        H1_Table_NormalizedAUPRC(row_id) = ...
            H1_NormalizedAUPRC(method_id, rep);
        H1_Table_TopK(row_id) = H1_TopKRecovery(method_id, rep);
        H1_Table_Budget95Edges(row_id) = H1_Budget95Edges(method_id, rep);
        H1_Table_Budget95Factor(row_id) = H1_Budget95Factor(method_id, rep);
        H1_Table_Precision95(row_id) = H1_PrecisionAt95(method_id, rep);
        H1_Table_FP95(row_id) = H1_FalsePositive95(method_id, rep);
        H1_Table_Pearson(row_id) = H1_PearsonWeightCorr(method_id, rep);
        H1_Table_Spearman(row_id) = H1_SpearmanWeightCorr(method_id, rep);
    end
end

H1_PerRepeat = table( ...
    H1_Table_Repeat, H1_Table_Method, ...
    H1_Table_TotalTargets, H1_Table_Priors, ...
    H1_Table_HeldoutTargets, H1_Table_Evaluated, H1_Table_Prevalence, ...
    H1_Table_AUPRC, H1_Table_NormalizedAUPRC, H1_Table_TopK, ...
    H1_Table_Budget95Edges, H1_Table_Budget95Factor, ...
    H1_Table_Precision95, H1_Table_FP95, ...
    H1_Table_Pearson, H1_Table_Spearman, ...
    'VariableNames', { ...
        'Repeat', 'Method', 'TotalH1Targets', 'PriorCount', ...
        'HeldoutH1Targets', 'EvaluatedNonpriorFeatures', ...
        'PositiveFraction', 'AUPRC', 'NormalizedAUPRC', ...
        'TopKRecovery', 'Budget95Edges', 'Budget95Factor', ...
        'PrecisionAt95Recall', 'FalsePositiveBurden95', ...
        'PearsonWeightCorrelation', 'SpearmanWeightCorrelation'});

H1_Summary = make_method_summary_table( ...
    H1_method_names, ...
    H1_AUPRC, ...
    H1_NormalizedAUPRC, ...
    H1_TopKRecovery, ...
    H1_Budget95Factor, ...
    H1_PrecisionAt95, ...
    H1_FalsePositive95, ...
    H1_PearsonWeightCorr, ...
    H1_SpearmanWeightCorr);

H1_CandidateRecovery = make_candidate_table( ...
    H1_method_names, candidate_multipliers_H1, H1_candidate_recovery);

H1_MeanPRCurve = make_pr_curve_table( ...
    H1_method_names, recall_grid, H1_PR_grid);

% ---------------------- Signal correlation tables ------------------------

NsignalRows = Nrep * Nestimate * Nreference;
Corr_Repeat    = zeros(NsignalRows, 1);
Corr_Estimate  = strings(NsignalRows, 1);
Corr_Reference = strings(NsignalRows, 1);
Corr_Value     = zeros(NsignalRows, 1);

row_id = 0;
for rep = 1:Nrep
    for e = 1:Nestimate
        for r = 1:Nreference
            row_id = row_id + 1;
            Corr_Repeat(row_id) = rep;
            Corr_Estimate(row_id) = estimate_names{e};
            Corr_Reference(row_id) = reference_names{r};
            Corr_Value(row_id) = SignalCorrelation(rep, e, r);
        end
    end
end

SignalCorrelation_PerRepeat = table( ...
    Corr_Repeat, Corr_Estimate, Corr_Reference, Corr_Value, ...
    'VariableNames', { ...
        'Repeat', 'EstimatedSignal', 'GTReferenceSignal', 'Correlation'});

NsummaryRows = Nestimate * Nreference;
Summary_Estimate  = strings(NsummaryRows, 1);
Summary_Reference = strings(NsummaryRows, 1);
Summary_Mean      = zeros(NsummaryRows, 1);
Summary_SD        = zeros(NsummaryRows, 1);

row_id = 0;
for e = 1:Nestimate
    for r = 1:Nreference
        row_id = row_id + 1;
        values = squeeze(SignalCorrelation(:, e, r));
        Summary_Estimate(row_id) = estimate_names{e};
        Summary_Reference(row_id) = reference_names{r};
        Summary_Mean(row_id) = mean(values, 'omitnan');
        Summary_SD(row_id) = std(values, 0, 'omitnan');
    end
end

SignalCorrelation_Summary = table( ...
    Summary_Estimate, Summary_Reference, Summary_Mean, Summary_SD, ...
    'VariableNames', { ...
        'EstimatedSignal', 'GTReferenceSignal', ...
        'CorrelationMean', 'CorrelationSD'});

% ---------------------- Signal composition tables ------------------------

NcompositionRows = Nrep * Nestimate;
Comp_Repeat   = zeros(NcompositionRows, 1);
Comp_Estimate = strings(NcompositionRows, 1);
Comp_H1       = zeros(NcompositionRows, 1);
Comp_Hother   = zeros(NcompositionRows, 1);
Comp_Global   = zeros(NcompositionRows, 1);
Comp_BetaH1     = zeros(NcompositionRows, 1);
Comp_BetaHother = zeros(NcompositionRows, 1);
Comp_BetaGlobal = zeros(NcompositionRows, 1);
Comp_R2         = zeros(NcompositionRows, 1);

row_id = 0;
for rep = 1:Nrep
    for e = 1:Nestimate
        row_id = row_id + 1;
        Comp_Repeat(row_id) = rep;
        Comp_Estimate(row_id) = estimate_names{e};
        Comp_H1(row_id) = SignalProportion(rep, e, 1);
        Comp_Hother(row_id) = SignalProportion(rep, e, 2);
        Comp_Global(row_id) = SignalProportion(rep, e, 3);
        Comp_BetaH1(row_id) = SignalBeta(rep, e, 1);
        Comp_BetaHother(row_id) = SignalBeta(rep, e, 2);
        Comp_BetaGlobal(row_id) = SignalBeta(rep, e, 3);
        Comp_R2(row_id) = SignalR2(rep, e);
    end
end

SignalComposition_PerRepeat = table( ...
    Comp_Repeat, Comp_Estimate, ...
    Comp_H1, Comp_Hother, Comp_Global, ...
    Comp_BetaH1, Comp_BetaHother, Comp_BetaGlobal, Comp_R2, ...
    'VariableNames', { ...
        'Repeat', 'EstimatedSignal', ...
        'TargetH1Proportion', 'OtherHiddenProportion', ...
        'GlobalProportion', 'BetaTargetH1', ...
        'BetaOtherHidden', 'BetaGlobal', 'R2'});

CompSummary_Estimate = string(estimate_names);
CompSummary_H1_Mean = zeros(Nestimate, 1);
CompSummary_H1_SD = zeros(Nestimate, 1);
CompSummary_Hother_Mean = zeros(Nestimate, 1);
CompSummary_Hother_SD = zeros(Nestimate, 1);
CompSummary_Global_Mean = zeros(Nestimate, 1);
CompSummary_Global_SD = zeros(Nestimate, 1);
CompSummary_R2_Mean = zeros(Nestimate, 1);
CompSummary_R2_SD = zeros(Nestimate, 1);

for e = 1:Nestimate
    v1 = squeeze(SignalProportion(:, e, 1));
    v2 = squeeze(SignalProportion(:, e, 2));
    v3 = squeeze(SignalProportion(:, e, 3));
    vr = squeeze(SignalR2(:, e));

    CompSummary_H1_Mean(e) = mean(v1, 'omitnan');
    CompSummary_H1_SD(e) = std(v1, 0, 'omitnan');
    CompSummary_Hother_Mean(e) = mean(v2, 'omitnan');
    CompSummary_Hother_SD(e) = std(v2, 0, 'omitnan');
    CompSummary_Global_Mean(e) = mean(v3, 'omitnan');
    CompSummary_Global_SD(e) = std(v3, 0, 'omitnan');
    CompSummary_R2_Mean(e) = mean(vr, 'omitnan');
    CompSummary_R2_SD(e) = std(vr, 0, 'omitnan');
end

SignalComposition_Summary = table( ...
    CompSummary_Estimate, ...
    CompSummary_H1_Mean, CompSummary_H1_SD, ...
    CompSummary_Hother_Mean, CompSummary_Hother_SD, ...
    CompSummary_Global_Mean, CompSummary_Global_SD, ...
    CompSummary_R2_Mean, CompSummary_R2_SD, ...
    'VariableNames', { ...
        'EstimatedSignal', ...
        'TargetH1Proportion_Mean', 'TargetH1Proportion_SD', ...
        'OtherHiddenProportion_Mean', 'OtherHiddenProportion_SD', ...
        'GlobalProportion_Mean', 'GlobalProportion_SD', ...
        'R2_Mean', 'R2_SD'});

% ----------------------------- Diagnostics -------------------------------

Repeat = (1:Nrep)';
Diagnostics = table( ...
    Repeat, ...
    LV_RetainedRank, ...
    LV_OracleSelectedComponent, ...
    LV_OracleCorrelationH1, ...
    LV_PriorSelectedComponent, ...
    LV_PriorEnrichment, ...
    GT_HiddenBlockOffdiagRatio, ...
    GT_HiddenGlobalCouplingNorm, ...
    'VariableNames', { ...
        'Repeat', ...
        'LVRetainedCanonicalRank', ...
        'LVOracleSelectedComponent', ...
        'LVOracleCorrelationWithGTH1', ...
        'LVPriorSelectedComponent', ...
        'LVPriorEnrichment', ...
        'GTHiddenBlockOffdiagRatio', ...
        'GTHiddenGlobalCouplingNorm'});

%% ========================================================================
% Combined H1 recovery and signal-composition figure
% Fixed normalized positions preserve the three panel dimensions.
%% ========================================================================

mean_prop = squeeze(mean(SignalProportion, 1, 'omitnan'));

fig2 = figure( ...
    'Color', 'w', ...
    'Units', 'centimeters', ...
    'Position', [2, 2, figure2_width_cm, figure2_height_cm], ...
    'PaperUnits', 'centimeters', ...
    'PaperSize', [figure2_width_cm, figure2_height_cm], ...
    'PaperPosition', [0, 0, figure2_width_cm, figure2_height_cm], ...
    'PaperPositionMode', 'manual');

% Equal panel widths within the 178 mm figure.
panel_y = 0.205;
panel_h = 0.715;
panel_w = 0.265;
ax1 = axes('Parent', fig2, 'Units', 'normalized', ...
    'Position', [0.065, panel_y, panel_w, panel_h]);
ax2 = axes('Parent', fig2, 'Units', 'normalized', ...
    'Position', [0.385, panel_y, panel_w, panel_h]);
ax3 = axes('Parent', fig2, 'Units', 'normalized', ...
    'Position', [0.705, panel_y, panel_w, panel_h]);

% ---------------- Panel A: held-out H1 PR curves -------------------------
hold(ax1, 'on');
pr_handles = gobjects(Nmethod_H1, 1);

for method_id = 1:Nmethod_H1

    light_color = 0.72 * ones(1, 3) + ...
        0.28 * H1_method_colors(method_id, :);

    for rep = 1:Nrep
        plot(ax1, recall_grid, H1_PR_grid(:, method_id, rep), ...
            'LineWidth', 0.65, ...
            'Color', light_color, ...
            'HandleVisibility', 'off');
    end

    mean_curve = mean( ...
        squeeze(H1_PR_grid(:, method_id, :)), ...
        2, 'omitnan');

    pr_handles(method_id) = plot(ax1, recall_grid, mean_curve, ...
        'LineWidth', 1.8, ...
        'Color', H1_method_colors(method_id, :), ...
        'DisplayName', H1_method_names{method_id});
end

baseline_H1 = mean(H1_PositiveFraction, 'omitnan');
yline(ax1, baseline_H1, ':', ...
    'Color', [0.35, 0.35, 0.35], ...
    'LineWidth', 0.8, ...
    'HandleVisibility', 'off');

xlim(ax1, [0, 1]);
ylim(ax1, [0, 1.015]);
xticks(ax1, [0, 0.5, 1]);
yticks(ax1, 0:0.2:1);
xlabel(ax1, 'Recall');
ylabel(ax1, 'Precision');
title(ax1, '(A) Held-out H1 recovery', ...
    'FontWeight', 'bold', ...
    'FontSize', font_size_title);
grid(ax1, 'on');
box(ax1, 'off');
text(ax1, 0.54, baseline_H1 + 0.026, 'Random baseline', ...
    'FontName', font_name, ...
    'FontSize', font_size_annotation, ...
    'Color', [0.2, 0.2, 0.2], ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'bottom');
set(ax1, ...
    'FontName', font_name, ...
    'FontSize', font_size_axis, ...
    'LineWidth', 0.8, ...
    'Layer', 'top');
hold(ax1, 'off');

% ---------------- Panel B: candidate-list budgets -----------------------
hold(ax2, 'on');
budget_handles = gobjects(Nmethod_H1, 1);

for method_id = 1:Nmethod_H1
    recovery_values = squeeze( ...
        H1_candidate_recovery(:, method_id, :));

    mean_recovery = mean(recovery_values, 2, 'omitnan');
    sd_recovery = std(recovery_values, 0, 2, 'omitnan');

    budget_handles(method_id) = errorbar( ...
        ax2, candidate_multipliers_H1, ...
        mean_recovery, sd_recovery, ...
        'Color', H1_method_colors(method_id, :), ...
        'LineWidth', 1.4, ...
        'Marker', H1_method_markers{method_id}, ...
        'MarkerSize', 4.8, ...
        'MarkerFaceColor', 'w', ...
        'CapSize', 3, ...
        'DisplayName', H1_method_names{method_id});
end

xline(ax2, 1, '--', ...
    'Color', [0.35, 0.35, 0.35], ...
    'LineWidth', 0.8, ...
    'HandleVisibility', 'off');

yline(ax2, target_recall, ':', ...
    'Color', [0.35, 0.35, 0.35], ...
    'LineWidth', 0.8, ...
    'HandleVisibility', 'off');

set(ax2, 'XScale', 'log');
xlim(ax2, [0.42, 11.5]);
ylim(ax2, [0, 1.015]);
xticks(ax2, candidate_multipliers_H1);
xticklabels(ax2, string(candidate_multipliers_H1));
yticks(ax2, 0:0.2:1);
xlabel(ax2, 'Candidate-list size / K');
ylabel(ax2, 'Held-out H1 recovery');
text(ax2, 1.03, 0.35, 'K budget', ...
    'Rotation', 90, ...
    'FontName', font_name, ...
    'FontSize', font_size_annotation, ...
    'Color', [0.2, 0.2, 0.2], ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'bottom');
text(ax2, 3.0, target_recall - 0.018, '95% recovery', ...
    'FontName', font_name, ...
    'FontSize', font_size_annotation, ...
    'Color', [0.2, 0.2, 0.2], ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top');
title(ax2, '(B) Candidate-list recovery', ...
    'FontWeight', 'bold', ...
    'FontSize', font_size_title);
grid(ax2, 'on');
box(ax2, 'off');
set(ax2, ...
    'FontName', font_name, ...
    'FontSize', font_size_axis, ...
    'LineWidth', 0.8, ...
    'Layer', 'top');
hold(ax2, 'off');

% ---------------- Panel C: signal composition ---------------------------
hold(ax3, 'on');
composition_handles = bar(ax3, mean_prop, 'stacked', ...
    'BarWidth', 0.72);

composition_colors = [ ...
    0.3010, 0.7450, 0.9330; ...
    0.9290, 0.6940, 0.1250; ...
    0.5000, 0.5000, 0.5000];

for r = 1:Nreference
    composition_handles(r).FaceColor = composition_colors(r, :);
end

xticks(ax3, 1:Nestimate);
xticklabels(ax3, {'HR-GGM', 'LV-GGM'});
xtickangle(ax3, 0);
xlim(ax3, [0.45, Nestimate + 0.55]);
ylim(ax3, [0, 1]);
yticks(ax3, 0:0.2:1);
ylabel(ax3, 'Regression-based proportion');
title(ax3, '(C) H1 signal composition', ...
    'FontWeight', 'bold', ...
    'FontSize', font_size_title);
grid(ax3, 'on');
box(ax3, 'off');
set(ax3, ...
    'FontName', font_name, ...
    'FontSize', font_size_axis, ...
    'LineWidth', 0.8, ...
    'Layer', 'top');

for e = 1:Nestimate
    cumulative_bottom = 0;
    for r = 1:Nreference
        value = mean_prop(e, r);
        if isfinite(value) && value >= 0.055
            text(ax3, e, cumulative_bottom + value / 2, ...
                sprintf('%.2f', value), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'Color', 'k', ...
                'FontWeight', 'bold', ...
                'FontSize', font_size_annotation, ...
                'FontName', font_name);
        end
        cumulative_bottom = cumulative_bottom + value;
    end
end
hold(ax3, 'off');

% ---------------- Legend areas ------------------------------------------
% Separate axes place legends outside the data panels.
ax_method_legend = axes('Parent', fig2, 'Units', 'normalized', ...
    'Position', [0.07, 0.015, 0.58, 0.105], ...
    'Visible', 'off');
hold(ax_method_legend, 'on');
method_legend_handles = gobjects(Nmethod_H1, 1);
for method_id = 1:Nmethod_H1
    method_legend_handles(method_id) = plot(ax_method_legend, nan, nan, ...
        'LineWidth', 1.8, ...
        'Color', H1_method_colors(method_id, :), ...
        'Marker', H1_method_markers{method_id}, ...
        'MarkerSize', 4.5, ...
        'MarkerFaceColor', 'w');
end
lgd_methods = legend(ax_method_legend, method_legend_handles, ...
    {'HR-GGM', 'LV-GGM', 'Prior correlation'}, ...
    'Orientation', 'horizontal', ...
    'Location', 'north', ...
    'Box', 'off', ...
    'FontName', font_name, ...
    'FontSize', font_size_legend);
lgd_methods.ItemTokenSize = [14, 8];
hold(ax_method_legend, 'off');

ax_comp_legend = axes('Parent', fig2, 'Units', 'normalized', ...
    'Position', [0.64, 0.015, 0.33, 0.105], ...
    'Visible', 'off');
hold(ax_comp_legend, 'on');
component_legend_handles = gobjects(Nreference, 1);
for r = 1:Nreference
    component_legend_handles(r) = patch(ax_comp_legend, ...
        nan, nan, composition_colors(r, :), ...
        'EdgeColor', 'k');
end
lgd_comp = legend(ax_comp_legend, component_legend_handles, ...
    {'Target H1', 'Other hidden', 'Global'}, ...
    'Orientation', 'horizontal', ...
    'NumColumns', 3, ...
    'Location', 'north', ...
    'Box', 'off', ...
    'FontName', font_name, ...
    'FontSize', font_size_legend);
lgd_comp.ItemTokenSize = [12, 8];
hold(ax_comp_legend, 'off');

% Align both legends on the same horizontal row.
drawnow;
set(lgd_methods, 'Units', 'normalized');
set(lgd_comp, 'Units', 'normalized');
legend_y = 0.025;
legend_h = 0.055;
set(lgd_methods, 'Position', [0.145, legend_y, 0.475, legend_h]);
set(lgd_comp,    'Position', [0.685, legend_y, 0.285, legend_h]);

save_figure_pair(fig2, output_dir, ...
    'H1_observed_edge_recovery_and_signal_specificity_1x3', ...
    figure_resolution);
close(fig2);

%% ============================== Save Excel ===============================

excel_file = fullfile(output_dir, ...
    'H1_observed_edge_analysis.xlsx');

if exist(excel_file, 'file')
    delete(excel_file);
end

writetable(LL_PerRepeat, excel_file, 'Sheet', 'Observed_PerRepeat');
writetable(LL_Summary, excel_file, 'Sheet', 'Observed_Summary');
writetable(LL_CandidateRecovery, excel_file, ...
    'Sheet', 'Observed_CandidateRecovery');
writetable(LL_MeanPRCurve, excel_file, 'Sheet', 'Observed_MeanPRCurve');

writetable(H1_PerRepeat, excel_file, 'Sheet', 'H1_PerRepeat');
writetable(H1_Summary, excel_file, 'Sheet', 'H1_Summary');
writetable(H1_CandidateRecovery, excel_file, ...
    'Sheet', 'H1_CandidateRecovery');
writetable(H1_MeanPRCurve, excel_file, 'Sheet', 'H1_MeanPRCurve');

writetable(SignalCorrelation_PerRepeat, excel_file, ...
    'Sheet', 'SignalCorr_PerRepeat');
writetable(SignalCorrelation_Summary, excel_file, ...
    'Sheet', 'SignalCorr_Summary');

writetable(SignalComposition_PerRepeat, excel_file, ...
    'Sheet', 'SignalComp_PerRepeat');
writetable(SignalComposition_Summary, excel_file, ...
    'Sheet', 'SignalComp_Summary');

writetable(Diagnostics, excel_file, 'Sheet', 'Diagnostics');

%% ========================================================================
% Local functions
%% ========================================================================

function file_path = resolve_data_file(folder_path, file_prefix, repeat_id)

    base_name = sprintf('%s_rep%d', file_prefix, repeat_id);

    candidate_paths = { ...
        fullfile(folder_path, [base_name, '.txt']), ...
        fullfile(folder_path, [base_name, '.csv']), ...
        fullfile(folder_path, [base_name, '.dat']), ...
        fullfile(folder_path, base_name)};

    file_path = '';

    for i = 1:numel(candidate_paths)
        if exist(candidate_paths{i}, 'file')
            file_path = candidate_paths{i};
            return;
        end
    end

    possible_files = dir(fullfile(folder_path, [base_name, '.*']));
    possible_files = possible_files(~[possible_files.isdir]);

    if numel(possible_files) == 1
        file_path = fullfile( ...
            possible_files(1).folder, possible_files(1).name);
        return;
    elseif numel(possible_files) > 1
        names = string({possible_files.name});
        error([ ...
            'Multiple files match "%s" in folder:\n%s\n%s'], ...
            base_name, folder_path, strjoin(names, newline));
    end

    error('Cannot find "%s" in folder:\n%s', base_name, folder_path);
end

function A = read_numeric_array(file_path)

    try
        A = readmatrix(file_path);
    catch first_error
        try
            A = dlmread(file_path);
        catch
            rethrow(first_error);
        end
    end

    if isempty(A)
        error('The file is empty: %s', file_path);
    end

    valid_rows = ~all(isnan(A), 2);
    valid_cols = ~all(isnan(A), 1);
    A = A(valid_rows, valid_cols);

    if isempty(A)
        error('No numeric data were detected in: %s', file_path);
    end

    A = double(A);

    if any(~isfinite(A(:)))
        error('Non-finite values were detected in: %s', file_path);
    end
end

function X_local = prepare_local_feature_data( ...
    Data_raw, Nlocal, allowed_feature_counts, repeat_id)

    [nrow, ncol] = size(Data_raw);
    allowed_feature_counts = unique(allowed_feature_counts(:)');

    row_matches = ismember(nrow, allowed_feature_counts);
    col_matches = ismember(ncol, allowed_feature_counts);

    if row_matches && ~col_matches
        % Features x samples
        X_local = Data_raw(1:Nlocal, :);

    elseif col_matches && ~row_matches
        % Samples x features
        X_local = Data_raw(:, 1:Nlocal)';

    elseif row_matches && col_matches
        error([ ...
            'Repeat %d observed-data size %d x %d is ambiguous because ', ...
            'both dimensions match an allowed feature count.'], ...
            repeat_id, nrow, ncol);

    elseif nrow >= Nlocal && nrow <= Nlocal + 20 && ncol > nrow
        % Accept feature dimensions near the expected matrix size.
        X_local = Data_raw(1:Nlocal, :);

    elseif ncol >= Nlocal && ncol <= Nlocal + 20 && nrow > ncol
        X_local = Data_raw(:, 1:Nlocal)';

    else
        error([ ...
            'Repeat %d observed-data size is %d x %d. One dimension must ', ...
            'equal a supported feature count (%s), or be between %d and ', ...
            '%d while the other dimension is the sample dimension.'], ...
            repeat_id, nrow, ncol, mat2str(allowed_feature_counts), ...
            Nlocal, Nlocal + 20);
    end

    if size(X_local, 1) ~= Nlocal
        error('Repeat %d local data must contain %d feature rows.', ...
            repeat_id, Nlocal);
    end

    if size(X_local, 2) < 3
        error('Repeat %d contains fewer than three samples.', repeat_id);
    end

    feature_sd = std(X_local, 0, 2);
    bad_features = find(~isfinite(feature_sd) | feature_sd <= eps);

    if ~isempty(bad_features)
        error([ ...
            'Repeat %d contains %d constant or invalid local features. ', ...
            'First affected indices: %s'], ...
            repeat_id, numel(bad_features), ...
            mat2str(bad_features(1:min(10, numel(bad_features)))'));
    end

end

function scores_full = calculate_prior_correlation_scores( ...
    X_local, prior_idx_matlab)

    Nlocal = size(X_local, 1);

    if any(prior_idx_matlab < 1) || any(prior_idx_matlab > Nlocal)
        error('Prior indices are outside the local-feature range.');
    end
    % Compute the Pearson correlation matrix across observed variables.
    R = corr(X_local', ...
        'Type', 'Pearson', ...
        'Rows', 'complete');

    if ~isequal(size(R), [Nlocal, Nlocal])
        error('Unexpected correlation-matrix size.');
    end

    if any(~isfinite(R(:)))
        error([ ...
            'Non-finite Pearson correlations were detected. Check for ', ...
            'constant features or invalid observed data.']);
    end

    % Score each feature by its mean absolute correlation with the prior set.
    scores_full = mean(abs(R(:, prior_idx_matlab)), 2);

    % Set prior-feature scores to NaN because they are excluded from evaluation.
    scores_full(prior_idx_matlab) = NaN;
end

function idx_cpp = read_cpp_index_vector(file_path)

    raw = read_numeric_array(file_path);
    idx_cpp = raw(:);

    if any(abs(idx_cpp - round(idx_cpp)) > 1e-10)
        error('The C++ index file contains non-integer values: %s', ...
            file_path);
    end

    idx_cpp = round(idx_cpp);
end

function assert_matrix_size(A, expected_rows, expected_cols, matrix_name)

    if ~isequal(size(A), [expected_rows, expected_cols])
        error([ ...
            '%s must be %d x %d, but detected size is %d x %d.'], ...
            matrix_name, expected_rows, expected_cols, ...
            size(A, 1), size(A, 2));
    end
end

function A = symmetrize_matrix(A)

    if size(A, 1) ~= size(A, 2)
        error('symmetrize_matrix requires a square matrix.');
    end

    A = (A + A') / 2;
end

function L = latent_component_signal(Theta, idxL, idxZ)

    theta_lz = Theta(idxL, idxZ);
    theta_zl = Theta(idxZ, idxL);
    theta_zz = Theta(idxZ, idxZ);

    if ~isscalar(theta_zz) || ...
            ~isfinite(theta_zz) || abs(theta_zz) <= eps
        error('Invalid latent diagonal in latent_component_signal.');
    end

    L = theta_lz * (theta_zz \ theta_zl);
    L = symmetrize_matrix(L);
end

function Prepared = prepare_canonical_components( ...
    input_array, ...
    Nlocal, ...
    relative_tolerance, ...
    absolute_tolerance, ...
    symmetry_tolerance, ...
    auto_orient, ...
    method_name, ...
    repeat_id)

    A = input_array;

    % Accept factor matrices with observed variables stored by columns.
    if size(A, 1) ~= Nlocal && size(A, 2) == Nlocal
        A = A';

    end

    if size(A, 1) ~= Nlocal
        error([ ...
            '%s repeat %d latent input must contain %d local-feature ', ...
            'rows after optional transposition. Detected %d x %d.'], ...
            method_name, repeat_id, Nlocal, size(A, 1), size(A, 2));
    end

    is_square = size(A, 2) == Nlocal;

    if is_square
        asymmetry_ratio = norm(A - A', 'fro') / max(norm(A, 'fro'), eps);
    else
        asymmetry_ratio = NaN;
    end

    if is_square && asymmetry_ratio <= symmetry_tolerance
        input_type = "symmetric latent matrix";
        L = symmetrize_matrix(A);
    else
        input_type = "factor matrix reconstructed as B*B'";
        L = A * A';
        L = symmetrize_matrix(L);
    end

    [U, eigenvalues] = eig(L, 'vector');
    U = real(U);
    eigenvalues = real(eigenvalues);

    positive_mass = sum(eigenvalues(eigenvalues > 0));
    negative_mass = sum(abs(eigenvalues(eigenvalues < 0)));

    if auto_orient && negative_mass > positive_mass
        L = -L;

        [U, eigenvalues] = eig(L, 'vector');
        U = real(U);
        eigenvalues = real(eigenvalues);
    end

    [eigenvalues, order] = sort(eigenvalues, 'descend');
    U = U(:, order);

    maximum_positive = max(eigenvalues);

    if ~isfinite(maximum_positive) || maximum_positive <= 0
        error([ ...
            '%s repeat %d contains no positive canonical ', ...
            'latent eigenvalue.'], ...
            method_name, repeat_id);
    end

    threshold = max( ...
        absolute_tolerance, ...
        relative_tolerance * maximum_positive);

    keep = eigenvalues > threshold;
    eigenvalues = eigenvalues(keep);
    U = U(:, keep);

    n_components = numel(eigenvalues);

    if n_components < 1
        error('%s repeat %d has no retained canonical component.', ...
            method_name, repeat_id);
    end

    components = cell(n_components, 1);
    feature_scores = cell(n_components, 1);

    for k = 1:n_components
        components{k} = symmetrize_matrix( ...
            eigenvalues(k) * U(:, k) * U(:, k)');

        feature_scores{k} = ...
            sqrt(eigenvalues(k)) * abs(U(:, k));
    end

    Prepared = struct();
    Prepared.components = components;
    Prepared.feature_scores = feature_scores;
    Prepared.eigenvalues = eigenvalues(:);
    Prepared.U = U;
    Prepared.n_components = n_components;
    Prepared.input_type = input_type;
end

function Selection = select_oracle_component(Prepared, S_GT_H1)

    correlation_H1 = nan(Prepared.n_components, 1);

    for k = 1:Prepared.n_components
        correlation_H1(k) = upper_triangle_correlation( ...
            Prepared.components{k}, S_GT_H1);
    end

    correlation_H1(~isfinite(correlation_H1)) = -Inf;
    [best_correlation, selected_rank] = max(correlation_H1);

    if ~isfinite(best_correlation)
        error('No finite oracle H1 correlation was found.');
    end

    Selection = struct();
    Selection.selected_rank = selected_rank;
    Selection.selected_component = Prepared.components{selected_rank};
    Selection.selected_feature_score = ...
        Prepared.feature_scores{selected_rank};
    Selection.selected_correlation_H1 = best_correlation;
end

function Selection = select_component_by_prior(Prepared, prior_idx_matlab)

    prior_enrichment = nan(Prepared.n_components, 1);

    for k = 1:Prepared.n_components
        score = Prepared.feature_scores{k};
        score_z = zscore_safe(score);
        prior_enrichment(k) = ...
            mean(score_z(prior_idx_matlab), 'omitnan');
    end

    prior_enrichment(~isfinite(prior_enrichment)) = -Inf;
    [best_enrichment, selected_rank] = max(prior_enrichment);

    if ~isfinite(best_enrichment)
        error('No finite prior-enrichment score was found.');
    end

    Selection = struct();
    Selection.selected_rank = selected_rank;
    Selection.selected_component = Prepared.components{selected_rank};
    Selection.selected_feature_score = ...
        Prepared.feature_scores{selected_rank};
    Selection.selected_prior_enrichment = best_enrichment;
end

function r = upper_triangle_correlation(A, B)

    if ~isequal(size(A), size(B))
        error('Matrices have different sizes in upper_triangle_correlation.');
    end

    n = size(A, 1);
    upper_mask = triu(true(n), 1);

    x = A(upper_mask);
    y = B(upper_mask);

    r = safe_corr(x, y, 'Pearson');
end

function output = calculate_ranking_metrics( ...
    scores, labels, K, target_recall, candidate_multipliers, recall_grid)

    scores = scores(:);
    labels = logical(labels(:));

    if numel(scores) ~= numel(labels)
        error('Score and label vectors have different lengths.');
    end

    if K <= 0
        error('K must be positive.');
    end

    scores(~isfinite(scores)) = 0;
    N = numel(scores);

    [sorted_scores, order] = sort(scores, 'descend');
    sorted_labels = labels(order);

    % Evaluate exact score ties as groups.
    group_start = [1; find(diff(sorted_scores) ~= 0) + 1];
    group_end = [group_start(2:end) - 1; N];

    Ngroup = numel(group_start);
    group_size = zeros(Ngroup, 1);
    group_positive = zeros(Ngroup, 1);

    for g = 1:Ngroup
        current_indices = group_start(g):group_end(g);
        group_size(g) = numel(current_indices);
        group_positive(g) = sum(sorted_labels(current_indices));
    end

    cumulative_selected = cumsum(group_size);
    cumulative_positive = cumsum(group_positive);

    recall = cumulative_positive / K;
    precision = cumulative_positive ./ cumulative_selected;

    delta_recall = diff([0; recall]);
    average_precision = sum(delta_recall .* precision);

    expected_TP_at_K = expected_positive_at_budget( ...
        group_size, group_positive, K);
    topK_recovery = expected_TP_at_K / K;

    target_positive_count = ceil(target_recall * K);
    target_budget = expected_budget_for_positive_count( ...
        group_size, group_positive, target_positive_count);

    target_budget_factor = target_budget / K;
    precision_at_target = target_positive_count / target_budget;
    false_positive_burden = target_budget - target_positive_count;

    candidate_recovery = zeros(numel(candidate_multipliers), 1);

    for b = 1:numel(candidate_multipliers)
        candidate_budget = round(candidate_multipliers(b) * K);
        candidate_budget = max(1, candidate_budget);
        candidate_budget = min(N, candidate_budget);

        expected_TP = expected_positive_at_budget( ...
            group_size, group_positive, candidate_budget);
        candidate_recovery(b) = expected_TP / K;
    end

    recall_curve = [0; recall];
    precision_curve = [1; precision];
    precision_envelope = flipud(cummax(flipud(precision_curve)));

    precision_on_grid = zeros(numel(recall_grid), 1);

    for q = 1:numel(recall_grid)
        requested_recall = recall_grid(q);
        curve_index = find(recall_curve >= requested_recall, 1, 'first');

        if isempty(curve_index)
            precision_on_grid(q) = precision_curve(end);
        else
            precision_on_grid(q) = precision_envelope(curve_index);
        end
    end

    output = struct();
    output.AUPRC = average_precision;
    output.TopKRecovery = topK_recovery;
    output.TargetBudgetEdges = target_budget;
    output.TargetBudgetFactor = target_budget_factor;
    output.PrecisionAtTarget = precision_at_target;
    output.FalsePositiveBurden = false_positive_burden;
    output.CandidateRecovery = candidate_recovery;
    output.PrecisionOnRecallGrid = precision_on_grid;
end

function expected_positive = expected_positive_at_budget( ...
    group_size, group_positive, budget)

    cumulative_size = cumsum(group_size);
    cumulative_positive = cumsum(group_positive);

    total_size = cumulative_size(end);
    budget = min(max(budget, 0), total_size);

    if budget == 0
        expected_positive = 0;
        return;
    end

    group_id = find(cumulative_size >= budget, 1, 'first');

    if group_id == 1
        previous_size = 0;
        previous_positive = 0;
    else
        previous_size = cumulative_size(group_id - 1);
        previous_positive = cumulative_positive(group_id - 1);
    end

    selected_from_tie = budget - previous_size;
    positive_fraction_in_tie = ...
        group_positive(group_id) / group_size(group_id);

    expected_positive = previous_positive + ...
        selected_from_tie * positive_fraction_in_tie;
end

function expected_budget = expected_budget_for_positive_count( ...
    group_size, group_positive, target_positive)

    cumulative_size = cumsum(group_size);
    cumulative_positive = cumsum(group_positive);

    group_id = find(cumulative_positive >= target_positive, 1, 'first');

    if isempty(group_id)
        expected_budget = cumulative_size(end);
        return;
    end

    if group_id == 1
        previous_size = 0;
        previous_positive = 0;
    else
        previous_size = cumulative_size(group_id - 1);
        previous_positive = cumulative_positive(group_id - 1);
    end

    positive_needed = target_positive - previous_positive;
    positive_fraction_in_tie = ...
        group_positive(group_id) / group_size(group_id);

    if positive_fraction_in_tie <= 0
        expected_budget = cumulative_size(group_id);
    else
        expected_selected_from_tie = ...
            positive_needed / positive_fraction_in_tie;
        expected_selected_from_tie = min( ...
            expected_selected_from_tie, group_size(group_id));
        expected_budget = previous_size + expected_selected_from_tie;
    end
end

function value = safe_corr(x, y, corr_type)

    x = x(:);
    y = y(:);

    keep = isfinite(x) & isfinite(y);
    x = x(keep);
    y = y(keep);

    if numel(x) < 3 || std(x) < eps || std(y) < eps
        value = NaN;
        return;
    end

    value = corr(x, y, ...
        'Type', corr_type, ...
        'Rows', 'complete');
end

function result = standardized_component_regression(y, Xcomponents)

    y = y(:);
    Xcomponents = double(Xcomponents);

    keep = isfinite(y) & all(isfinite(Xcomponents), 2);
    y = y(keep);
    Xcomponents = Xcomponents(keep, :);

    yz = zscore_safe(y);

    Xz = zeros(size(Xcomponents));
    for j = 1:size(Xcomponents, 2)
        Xz(:, j) = zscore_safe(Xcomponents(:, j));
    end

    X = [ones(size(Xz, 1), 1), Xz];

    if rcond(X' * X) < 1e-12
        beta_full = pinv(X) * yz;
    else
        beta_full = X \ yz;
    end

    yhat = X * beta_full;

    denominator = sum((yz - mean(yz)).^2);
    if denominator < eps
        R2 = NaN;
    else
        R2 = 1 - sum((yz - yhat).^2) / denominator;
    end

    beta_components = beta_full(2:end);
    abs_beta = abs(beta_components);
    beta_sum = sum(abs_beta);

    if beta_sum < eps
        proportion = nan(size(abs_beta));
    else
        proportion = abs_beta / beta_sum;
    end

    result = struct();
    result.beta = beta_components(:);
    result.proportion = proportion(:);
    result.R2 = R2;
end

function z = zscore_safe(x)

    x = x(:);
    mu = mean(x, 'omitnan');
    sigma = std(x, 0, 'omitnan');

    if ~isfinite(sigma) || sigma < eps
        z = zeros(size(x));
    else
        z = (x - mu) / sigma;
    end
end

function T = make_method_summary_table( ...
    method_names, AUPRC, NormalizedAUPRC, TopKRecovery, ...
    Budget95Factor, PrecisionAt95, FalsePositive95, ...
    PearsonCorr, SpearmanCorr)

    Nmethod = numel(method_names);
    Method = string(method_names(:));

    AUPRC_Mean = zeros(Nmethod, 1);
    AUPRC_SD = zeros(Nmethod, 1);
    NormalizedAUPRC_Mean = zeros(Nmethod, 1);
    NormalizedAUPRC_SD = zeros(Nmethod, 1);
    TopKRecovery_Mean = zeros(Nmethod, 1);
    TopKRecovery_SD = zeros(Nmethod, 1);
    Budget95Factor_Mean = zeros(Nmethod, 1);
    Budget95Factor_SD = zeros(Nmethod, 1);
    PrecisionAt95Recall_Mean = zeros(Nmethod, 1);
    PrecisionAt95Recall_SD = zeros(Nmethod, 1);
    FalsePositiveBurden95_Mean = zeros(Nmethod, 1);
    FalsePositiveBurden95_SD = zeros(Nmethod, 1);

    for m = 1:Nmethod
        AUPRC_Mean(m) = mean(AUPRC(m, :), 'omitnan');
        AUPRC_SD(m) = std(AUPRC(m, :), 0, 'omitnan');
        NormalizedAUPRC_Mean(m) = ...
            mean(NormalizedAUPRC(m, :), 'omitnan');
        NormalizedAUPRC_SD(m) = ...
            std(NormalizedAUPRC(m, :), 0, 'omitnan');
        TopKRecovery_Mean(m) = ...
            mean(TopKRecovery(m, :), 'omitnan');
        TopKRecovery_SD(m) = ...
            std(TopKRecovery(m, :), 0, 'omitnan');
        Budget95Factor_Mean(m) = ...
            mean(Budget95Factor(m, :), 'omitnan');
        Budget95Factor_SD(m) = ...
            std(Budget95Factor(m, :), 0, 'omitnan');
        PrecisionAt95Recall_Mean(m) = ...
            mean(PrecisionAt95(m, :), 'omitnan');
        PrecisionAt95Recall_SD(m) = ...
            std(PrecisionAt95(m, :), 0, 'omitnan');
        FalsePositiveBurden95_Mean(m) = ...
            mean(FalsePositive95(m, :), 'omitnan');
        FalsePositiveBurden95_SD(m) = ...
            std(FalsePositive95(m, :), 0, 'omitnan');
    end

    T = table( ...
        Method, ...
        AUPRC_Mean, AUPRC_SD, ...
        NormalizedAUPRC_Mean, NormalizedAUPRC_SD, ...
        TopKRecovery_Mean, TopKRecovery_SD, ...
        Budget95Factor_Mean, Budget95Factor_SD, ...
        PrecisionAt95Recall_Mean, PrecisionAt95Recall_SD, ...
        FalsePositiveBurden95_Mean, FalsePositiveBurden95_SD);

    if ~isempty(PearsonCorr)
        PearsonWeightCorrelation_Mean = zeros(Nmethod, 1);
        PearsonWeightCorrelation_SD = zeros(Nmethod, 1);
        SpearmanWeightCorrelation_Mean = zeros(Nmethod, 1);
        SpearmanWeightCorrelation_SD = zeros(Nmethod, 1);

        for m = 1:Nmethod
            PearsonWeightCorrelation_Mean(m) = ...
                mean(PearsonCorr(m, :), 'omitnan');
            PearsonWeightCorrelation_SD(m) = ...
                std(PearsonCorr(m, :), 0, 'omitnan');
            SpearmanWeightCorrelation_Mean(m) = ...
                mean(SpearmanCorr(m, :), 'omitnan');
            SpearmanWeightCorrelation_SD(m) = ...
                std(SpearmanCorr(m, :), 0, 'omitnan');
        end

        T.PearsonWeightCorrelation_Mean = ...
            PearsonWeightCorrelation_Mean;
        T.PearsonWeightCorrelation_SD = ...
            PearsonWeightCorrelation_SD;
        T.SpearmanWeightCorrelation_Mean = ...
            SpearmanWeightCorrelation_Mean;
        T.SpearmanWeightCorrelation_SD = ...
            SpearmanWeightCorrelation_SD;
    end
end

function T = make_candidate_table( ...
    method_names, candidate_multipliers, recovery_array)

    Nmethod = numel(method_names);
    Nbudget = numel(candidate_multipliers);
    Nrep = size(recovery_array, 3);
    Nrow = Nmethod * Nbudget * Nrep;

    Repeat = zeros(Nrow, 1);
    Method = strings(Nrow, 1);
    CandidateBudgetTimesK = zeros(Nrow, 1);
    Recall = zeros(Nrow, 1);

    row_id = 0;
    for rep = 1:Nrep
        for method_id = 1:Nmethod
            for b = 1:Nbudget
                row_id = row_id + 1;
                Repeat(row_id) = rep;
                Method(row_id) = method_names{method_id};
                CandidateBudgetTimesK(row_id) = ...
                    candidate_multipliers(b);
                Recall(row_id) = ...
                    recovery_array(b, method_id, rep);
            end
        end
    end

    T = table(Repeat, Method, CandidateBudgetTimesK, Recall);
end

function T = make_pr_curve_table(method_names, recall_grid, PR_array)

    T = table(recall_grid, 'VariableNames', {'Recall'});

    for method_id = 1:numel(method_names)
        current_values = squeeze(PR_array(:, method_id, :));
        mean_curve = mean(current_values, 2, 'omitnan');
        sd_curve = std(current_values, 0, 2, 'omitnan');

        base_name = matlab.lang.makeValidName(method_names{method_id});
        T.([base_name, '_MeanPrecision']) = mean_curve;
        T.([base_name, '_SDPrecision']) = sd_curve;
    end
end

function save_figure_pair(fig_handle, output_dir, base_name, resolution)

    png_file = fullfile(output_dir, [base_name, '.png']);
    pdf_file = fullfile(output_dir, [base_name, '.pdf']);

    set(fig_handle, 'InvertHardcopy', 'off');
    drawnow;
    exportgraphics(fig_handle, png_file, ...
        'Resolution', resolution, ...
        'BackgroundColor', 'white');
    exportgraphics(fig_handle, pdf_file, ...
        'ContentType', 'vector', ...
        'BackgroundColor', 'white');
end
