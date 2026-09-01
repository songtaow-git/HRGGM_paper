%% HR-GGM parameter-selection simulation data generator
% Generates matched Gaussian simulation datasets for parameter-selection
% validation with 200 observed features, one hidden feature, one global
% feature, 10,000 samples, three repeats, and three noise levels.
%
% Each condition-specific Data folder contains:
%   Alpha_list.txt
%   Data_whole.bin
%   Idx_g_list.txt
%   Idx_h_list.txt
%   Lambda_list.txt
%   Theta_full.txt
%
% Data_whole.bin stores uint64 row and column counts followed by
% double-precision values for a local-feature-by-sample matrix.
% Overall_simulation_summary.csv is saved in the root output folder.

clear;
clc;
close all;

%% Configure output and parameter-list inputs
RootOutputFolder = "ParameterValidation_200L_1H_1G_10000_rep3";
reset_output_folder = true;

Alpha_list_file  = "Alpha_list.txt";
Lambda_list_file = "Lambda_list.txt";

if ~exist(Alpha_list_file, 'file')
    error('Cannot find Alpha_list_file: %s', Alpha_list_file);
end

if ~exist(Lambda_list_file, 'file')
    error('Cannot find Lambda_list_file: %s', Lambda_list_file);
end

if reset_output_folder && exist(RootOutputFolder, 'dir')
    rmdir(RootOutputFolder, 's');
end

if ~exist(RootOutputFolder, 'dir')
    mkdir(RootOutputFolder);
end

%% Define simulation dimensions
n_nodes_per_group = 5;
n_groups          = 40;
n_genes           = n_nodes_per_group * n_groups;  % 200 local features

n_cells = 10000;

n_hidden = 1;
n_global = 1;
n_repeats = 3;

if n_hidden ~= 1 || n_global ~= 1
    error('This script requires exactly one hidden and one global feature.');
end

p = n_genes + n_hidden + n_global;

idx_genes  = 1:n_genes;
idx_hidden = n_genes + 1;  % 201 in MATLAB
idx_global = n_genes + 2;  % 202 in MATLAB

%% Define network structure and noise conditions
% The LL network contains 200 within-group ring edges and 50
% between-group edges. Eight complete groups provide 40 LH edges.
n_between_LL_edges = 50;

n_LH_edges = 40;

if mod(n_LH_edges, n_nodes_per_group) ~= 0
    error('n_LH_edges must be divisible by n_nodes_per_group.');
end

n_hidden_target_groups = n_LH_edges / n_nodes_per_group;

n_hidden_idx_per_group = 2;

noise_variance_fraction_grid = [0, 0.20, 0.50];

if any(noise_variance_fraction_grid < 0) || ...
        any(noise_variance_fraction_grid > 1)
    error('Noise variance fractions must lie between 0 and 1.');
end

%% Define precision-matrix parameters
local_unit_diag  = 2;
local_unit_w_min = 0.4;
local_unit_w_max = 0.8;
local_unit_sign  = -1;

between_w_min = local_unit_w_min;
between_w_max = local_unit_w_max;
between_sign  = local_unit_sign;

hidden_diag  = 4;
hidden_w_min = 0.45;
hidden_w_max = 0.75;

global_diag  = 8;
global_shift = 0.1;
global_span  = 0.5;

spd_margin = 0.1;

standardize_cov_diag = true;
numerical_zero_threshold = 1e-10;

%% Precompute group labels and eligible between-group pairs
group_id = repelem((1:n_groups)', n_nodes_per_group);

all_local_gene_pairs = nchoosek(1:n_genes, 2);
is_between_pair = ...
    group_id(all_local_gene_pairs(:,1)) ~= ...
    group_id(all_local_gene_pairs(:,2));

between_gene_pairs = all_local_gene_pairs(is_between_pair, :);
n_possible_between_edges = size(between_gene_pairs, 1);

if n_between_LL_edges > n_possible_between_edges
    error('Requested between-group LL edges exceed possible pairs.');
end

if n_hidden_target_groups > n_groups
    error('The number of hidden-connected groups exceeds n_groups.');
end

if n_hidden_idx_per_group < 1 || ...
        n_hidden_idx_per_group > n_nodes_per_group
    error('Invalid n_hidden_idx_per_group.');
end

%% Generate repeat-specific structures shared across noise levels
% Edge identities, edge weights, and hidden-target groups are fixed
% within each repeat and reused for all three noise conditions.
between_edge_order_by_rep  = cell(n_repeats, 1);
between_edge_weight_by_rep = cell(n_repeats, 1);
LH_group_order_by_rep      = cell(n_repeats, 1);
LH_group_weight_by_rep     = cell(n_repeats, 1);

for rep = 1:n_repeats
    rng(900000 + rep, 'twister');

    between_edge_order_by_rep{rep} = ...
        randperm(n_possible_between_edges);

    between_edge_weight_by_rep{rep} = ...
        between_sign * ...
        (between_w_min + ...
        (between_w_max - between_w_min) * ...
        rand(n_between_LL_edges, 1));

    rng(910000 + rep, 'twister');

    LH_group_order_by_rep{rep} = randperm(n_groups);

    LH_group_weight_by_rep{rep} = ...
        -(hidden_w_min + ...
        (hidden_w_max - hidden_w_min) * ...
        rand(n_groups, n_nodes_per_group));
end

fprintf('============================================================\n');
fprintf('HR-GGM parameter-validation simulation\n');
fprintf('Local features: %d\n', n_genes);
fprintf('Hidden features: %d\n', n_hidden);
fprintf('Global features: %d\n', n_global);
fprintf('Training samples: %d\n', n_cells);
fprintf('Repeats: %d\n', n_repeats);
fprintf('Within-group LL edges: %d\n', n_groups * n_nodes_per_group);
fprintf('Between-group LL edges: %d\n', n_between_LL_edges);
fprintf('Total LL edges: %d\n', ...
    n_groups * n_nodes_per_group + n_between_LL_edges);
fprintf('True LH edges: %d\n', n_LH_edges);
fprintf('Prior features: %d\n', ...
    n_hidden_target_groups * n_hidden_idx_per_group);
fprintf('Held-out LH targets: %d\n', ...
    n_LH_edges - n_hidden_target_groups * n_hidden_idx_per_group);
fprintf('Noise levels: '); fprintf('%.0f%% ', ...
    100 * noise_variance_fraction_grid); fprintf('\n');
fprintf('============================================================\n');

overall_summary = table();

%% Generate matched simulation datasets
for rep = 1:n_repeats

    fprintf('\n############################################################\n');
    fprintf('Generating repeat %d of %d\n', rep, n_repeats);
    fprintf('############################################################\n');

    % Construct the local-local precision block from 40 five-node rings.
    rng(100000 + rep, 'twister');

    Theta_LL_unshifted = zeros(n_genes, n_genes);

    within_edge_list = zeros(n_groups * n_nodes_per_group, 3);
    within_edge_counter = 0;

    for g = 1:n_groups

        group_idx = ...
            ((g - 1) * n_nodes_per_group + 1):...
            (g * n_nodes_per_group);

        Unit_precision_LL = local_unit_diag * eye(n_nodes_per_group);

        for j = 1:(n_nodes_per_group - 1)

            edge_weight = local_unit_sign * ...
                (local_unit_w_min + ...
                (local_unit_w_max - local_unit_w_min) * rand());

            Unit_precision_LL(j, j + 1) = edge_weight;
            Unit_precision_LL(j + 1, j) = edge_weight;

            within_edge_counter = within_edge_counter + 1;
            within_edge_list(within_edge_counter, :) = ...
                [group_idx(j), group_idx(j + 1), edge_weight];
        end

        edge_weight = local_unit_sign * ...
            (local_unit_w_min + ...
            (local_unit_w_max - local_unit_w_min) * rand());

        Unit_precision_LL(1, n_nodes_per_group) = edge_weight;
        Unit_precision_LL(n_nodes_per_group, 1) = edge_weight;

        within_edge_counter = within_edge_counter + 1;
        within_edge_list(within_edge_counter, :) = ...
            [group_idx(1), group_idx(end), edge_weight];

        Theta_LL_unshifted(group_idx, group_idx) = Unit_precision_LL;
    end

    within_edge_list = within_edge_list(1:within_edge_counter, :);

    % Add the repeat-specific between-group LL edges.
    selected_between_ids = ...
        between_edge_order_by_rep{rep}(1:n_between_LL_edges);

    selected_between_pairs = ...
        between_gene_pairs(selected_between_ids, :);

    selected_between_weights = ...
        between_edge_weight_by_rep{rep};

    between_edge_list = zeros(n_between_LL_edges, 3);

    for e = 1:n_between_LL_edges

        gene_i = selected_between_pairs(e, 1);
        gene_j = selected_between_pairs(e, 2);
        edge_weight = selected_between_weights(e);

        Theta_LL_unshifted(gene_i, gene_j) = edge_weight;
        Theta_LL_unshifted(gene_j, gene_i) = edge_weight;

        between_edge_list(e, :) = [gene_i, gene_j, edge_weight];
    end

    % Assign hidden-local edges to eight complete feature groups.
    hidden_target_groups = sort(...
        LH_group_order_by_rep{rep}(1:n_hidden_target_groups));

    Theta_LH_unshifted = zeros(n_genes, 1);
    hidden_targets_1based = zeros(n_LH_edges, 1);

    target_counter = 0;

    for k = 1:n_hidden_target_groups

        g = hidden_target_groups(k);

        group_idx = ...
            ((g - 1) * n_nodes_per_group + 1):...
            (g * n_nodes_per_group);

        hidden_group_weights = ...
            LH_group_weight_by_rep{rep}(g, :).';

        Theta_LH_unshifted(group_idx, 1) = hidden_group_weights;

        hidden_targets_1based(...
            target_counter + (1:n_nodes_per_group)) = group_idx(:);

        target_counter = target_counter + n_nodes_per_group;
    end

    hidden_targets_1based = sort(hidden_targets_1based);

    % Assign a repeated global-loading pattern to all local groups.
    global_repeat_weights = ...
        (rand(1, n_nodes_per_group) - 0.5) * global_span;

    for j = 1:n_nodes_per_group
        if global_repeat_weights(j) < 0
            global_repeat_weights(j) = ...
                global_repeat_weights(j) - global_shift;
        elseif global_repeat_weights(j) > 0
            global_repeat_weights(j) = ...
                global_repeat_weights(j) + global_shift;
        else
            global_repeat_weights(j) = global_shift;
        end
    end

    Theta_LG_unshifted = zeros(n_genes, 1);

    for g = 1:n_groups
        group_idx = ...
            ((g - 1) * n_nodes_per_group + 1):...
            (g * n_nodes_per_group);

        Theta_LG_unshifted(group_idx, 1) = global_repeat_weights.';
    end

    global_targets_1based = find(Theta_LG_unshifted ~= 0);

    % Assemble the complete precision matrix.
    Theta_unshifted = zeros(p, p);

    Theta_unshifted(idx_genes, idx_genes) = Theta_LL_unshifted;

    Theta_unshifted(idx_genes, idx_hidden) = Theta_LH_unshifted;
    Theta_unshifted(idx_hidden, idx_genes) = Theta_LH_unshifted.';

    Theta_unshifted(idx_genes, idx_global) = Theta_LG_unshifted;
    Theta_unshifted(idx_global, idx_genes) = Theta_LG_unshifted.';

    Theta_unshifted(idx_hidden, idx_hidden) = hidden_diag;
    Theta_unshifted(idx_global, idx_global) = global_diag;

    Theta_unshifted(idx_hidden, idx_global) = 0;
    Theta_unshifted(idx_global, idx_hidden) = 0;

    Theta_unshifted = (Theta_unshifted + Theta_unshifted.') / 2;

    % Apply the minimum diagonal shift required by the SPD margin.
    min_eig_unshifted = min(eig(Theta_unshifted));
    diagonal_shift = max(0, spd_margin - min_eig_unshifted);

    Theta_shifted = ...
        Theta_unshifted + diagonal_shift * eye(p);

    Theta_shifted = (Theta_shifted + Theta_shifted.') / 2;

    min_eig_after_shift = min(eig(Theta_shifted));
    [~, chol_flag_shifted] = chol(Theta_shifted);

    if chol_flag_shifted ~= 0
        error('Theta_shifted is not SPD in repeat %d.', rep);
    end

    fprintf('Repeat %d minimum eigenvalue before shift: %.10f\n', ...
        rep, min_eig_unshifted);
    fprintf('Repeat %d diagonal shift: %.10f\n', ...
        rep, diagonal_shift);
    fprintf('Repeat %d minimum eigenvalue after shift: %.10f\n', ...
        rep, min_eig_after_shift);

    % Convert the precision matrix to covariance form and standardize
    % marginal variances when requested.
    Sigma_raw = Theta_shifted \ eye(p);
    Sigma_raw = (Sigma_raw + Sigma_raw.') / 2;

    if standardize_cov_diag

        sd_vec = sqrt(diag(Sigma_raw));

        if any(~isfinite(sd_vec)) || any(sd_vec <= 0)
            error('Invalid covariance diagonal in repeat %d.', rep);
        end

        Sigma = Sigma_raw ./ (sd_vec * sd_vec.');
        Sigma = (Sigma + Sigma.') / 2;

    else
        Sigma = Sigma_raw;
    end

    Theta = Sigma \ eye(p);
    Theta = (Theta + Theta.') / 2;
    Theta(abs(Theta) < numerical_zero_threshold) = 0;

    [~, chol_flag_theta] = chol(Theta);

    if chol_flag_theta ~= 0
        error('Final standardized Theta is not SPD in repeat %d.', rep);
    end

    min_eig_final_theta = min(eig(Theta));
    diag_sigma_final = diag(Sigma);

    Theta_LL = Theta(idx_genes, idx_genes);

    % Verify the realized LL support and edge counts.
    GT_LL_mask = abs(Theta_LL) > numerical_zero_threshold;
    GT_LL_mask(1:n_genes+1:end) = false;

    same_group_matrix = bsxfun(@eq, group_id, group_id.');

    GT_LL_within_mask  = GT_LL_mask & same_group_matrix;
    GT_LL_between_mask = GT_LL_mask & ~same_group_matrix;

    n_LL_within_actual  = nnz(triu(GT_LL_within_mask, 1));
    n_LL_between_actual = nnz(triu(GT_LL_between_mask, 1));
    n_LL_total_actual   = nnz(triu(GT_LL_mask, 1));

    if n_LL_within_actual ~= 200
        error('Expected 200 within-group LL edges; found %d.', ...
            n_LL_within_actual);
    end

    if n_LL_between_actual ~= 50
        error('Expected 50 between-group LL edges; found %d.', ...
            n_LL_between_actual);
    end

    if n_LL_total_actual ~= 250
        error('Expected 250 total LL edges; found %d.', ...
            n_LL_total_actual);
    end

    % Select two prior features from each hidden-connected group and
    % assign all remaining local features to the global index list.
    rng(120000 + rep, 'twister');

    prior_idx_1based = zeros(...
        n_hidden_target_groups * n_hidden_idx_per_group, 1);

    prior_counter = 0;

    for k = 1:n_hidden_target_groups

        g = hidden_target_groups(k);

        group_idx = ...
            ((g - 1) * n_nodes_per_group + 1):...
            (g * n_nodes_per_group);

        selected_local_positions = ...
            randperm(n_nodes_per_group, n_hidden_idx_per_group);

        selected_prior_idx = group_idx(selected_local_positions);

        prior_idx_1based(...
            prior_counter + (1:n_hidden_idx_per_group)) = ...
            selected_prior_idx(:);

        prior_counter = prior_counter + n_hidden_idx_per_group;
    end

    prior_idx_1based = sort(prior_idx_1based);

    heldout_hidden_targets_1based = ...
        setdiff(hidden_targets_1based, prior_idx_1based);

    all_local_idx_1based = (1:n_genes).';
    global_idx_1based = ...
        setdiff(all_local_idx_1based, prior_idx_1based);

    Idx_h_list = prior_idx_1based - 1;
    Idx_g_list = global_idx_1based - 1;

    hidden_targets_cpp = hidden_targets_1based - 1;
    hidden_heldout_cpp = heldout_hidden_targets_1based - 1;
    global_targets_cpp = global_targets_1based - 1;

    if numel(Idx_h_list) ~= 16
        error('Expected 16 prior indices; found %d.', numel(Idx_h_list));
    end

    if numel(hidden_heldout_cpp) ~= 24
        error('Expected 24 held-out LH targets; found %d.', ...
            numel(hidden_heldout_cpp));
    end

    if numel(Idx_h_list) + numel(Idx_g_list) ~= n_genes
        error('Idx_h_list and Idx_g_list do not cover all local features.');
    end

    if ~isempty(intersect(Idx_h_list, Idx_g_list))
        error('Idx_h_list and Idx_g_list overlap.');
    end

    % Compute LL partial-correlation diagnostics for the root summary.
    diag_LL = sqrt(diag(Theta_LL));
    PartialCorr_LL = -Theta_LL ./ (diag_LL * diag_LL.');

    within_upper_mask  = triu(GT_LL_within_mask, 1);
    between_upper_mask = triu(GT_LL_between_mask, 1);

    abs_pc_within  = abs(PartialCorr_LL(within_upper_mask));
    abs_pc_between = abs(PartialCorr_LL(between_upper_mask));

    % Generate clean local observations and a matched Gaussian noise matrix.
    [L_sigma, chol_flag_sigma] = chol(Sigma, 'lower');

    if chol_flag_sigma ~= 0
        error('Sigma is not SPD in repeat %d.', rep);
    end

    rng(130000 + rep, 'twister');

    Z_clean = randn(n_cells, p) * L_sigma.';

    Z_signal = standardize_columns(...
        Z_clean(:, idx_genes));

    clear Z_clean;

    rng(140000 + rep, 'twister');

    Noise_z = standardize_columns(...
        randn(n_cells, n_genes));

    % Mix the shared signal and noise matrices at each variance fraction.
    for iNoise = 1:numel(noise_variance_fraction_grid)

        noise_fraction = noise_variance_fraction_grid(iNoise);
        noise_percent = round(100 * noise_fraction);

        signal_weight = sqrt(1 - noise_fraction);
        noise_weight  = sqrt(noise_fraction);

        Z_observed = ...
            signal_weight * Z_signal + ...
            noise_weight  * Noise_z;

        Z_observed = ...
            standardize_columns(Z_observed);

        condition_name = sprintf(...
            'Data_200L_LL250_LH040_noise_%03dpct_rep_%02d', ...
            noise_percent, rep);

        OutputFolder = ...
            fullfile(RootOutputFolder, condition_name, 'Data');

        if ~exist(OutputFolder, 'dir')
            mkdir(OutputFolder);
        end

        fprintf('\n------------------------------------------------------------\n');
        fprintf('Saving condition: %s\n', condition_name);
        fprintf('Noise variance fraction: %.2f\n', noise_fraction);
        fprintf('Output folder: %s\n', OutputFolder);
        fprintf('------------------------------------------------------------\n');

        % Write the condition-specific ground truth, index lists, and
        % parameter grids to the Data folder.
        writematrix(Theta, ...
            fullfile(OutputFolder, "Theta_full.txt"), ...
            'Delimiter', 'tab');

        write_matrix_or_empty(...
            fullfile(OutputFolder, "Idx_h_list.txt"), ...
            Idx_h_list);

        write_matrix_or_empty(...
            fullfile(OutputFolder, "Idx_g_list.txt"), ...
            Idx_g_list);

        copyfile(Alpha_list_file, ...
            fullfile(OutputFolder, "Alpha_list.txt"));

        copyfile(Lambda_list_file, ...
            fullfile(OutputFolder, "Lambda_list.txt"));

        % Store the observed matrix in the binary layout used by C++.
        Matrix = Z_observed.';  % 200 local features x 10,000 samples

        write_binary_matrix(...
            fullfile(OutputFolder, "Data_whole.bin"), ...
            Matrix);

        % Assemble condition-level diagnostics for the root summary.
        if isempty(abs_pc_within)
            pc_within_min = NaN;
            pc_within_median = NaN;
            pc_within_max = NaN;
        else
            pc_within_min = min(abs_pc_within);
            pc_within_median = median(abs_pc_within);
            pc_within_max = max(abs_pc_within);
        end

        if isempty(abs_pc_between)
            pc_between_min = NaN;
            pc_between_median = NaN;
            pc_between_max = NaN;
        else
            pc_between_min = min(abs_pc_between);
            pc_between_median = median(abs_pc_between);
            pc_between_max = max(abs_pc_between);
        end

        fprintf('Data_whole.bin matrix: %d x %d\n', ...
            size(Matrix,1), size(Matrix,2));
        fprintf('LL edges: %d within + %d between = %d total\n', ...
            n_LL_within_actual, ...
            n_LL_between_actual, ...
            n_LL_total_actual);
        fprintf('LH targets: %d total, %d prior, %d held out\n', ...
            n_LH_edges, numel(Idx_h_list), numel(hidden_heldout_cpp));
        fprintf('LG targets: %d\n', numel(global_targets_cpp));

        clear Z_observed Matrix
    end

    clear Theta_unshifted Theta_shifted Theta Sigma Sigma_raw
    clear Theta_LL Theta_LL_unshifted
    clear Z_signal Noise_z
end

%% Save the aggregate simulation summary

fprintf('\n============================================================\n');
fprintf('All simulations finished successfully.\n');
fprintf('Root output folder: %s\n', RootOutputFolder);
fprintf('Generated conditions: %d\n', height(overall_summary));
fprintf('Expected conditions: %d repeats x %d noise levels = %d\n', ...
    n_repeats, numel(noise_variance_fraction_grid), ...
    n_repeats * numel(noise_variance_fraction_grid));
fprintf('============================================================\n');

%% Local functions

% Center and scale each column using the sample standard deviation.
function Xz = standardize_columns(X)

    mu = mean(X, 1);
    sigma = std(X, 0, 1);

    sigma(~isfinite(sigma) | sigma == 0) = 1;

    Xz = (X - mu) ./ sigma;

    if any(~isfinite(Xz), 'all')
        error('Standardized matrix contains NaN or Inf.');
    end
end

% Write a feature-by-sample matrix in the C++ RowMat binary layout.
function write_binary_matrix(bin_file, Matrix)

    fid = fopen(bin_file, 'wb');

    if fid == -1
        error('Cannot open binary output file: %s', bin_file);
    end

    cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fwrite(fid, uint64(size(Matrix, 1)), 'uint64');
    fwrite(fid, uint64(size(Matrix, 2)), 'uint64');

    fwrite(fid, Matrix.', 'double');
end

% Write a tab-delimited vector or create an empty file.
function write_matrix_or_empty(filename, X)

    if isempty(X)

        fid = fopen(filename, 'w');

        if fid == -1
            error('Cannot create empty file: %s', filename);
        end

        fclose(fid);

    else
        writematrix(X, filename, 'Delimiter', 'tab');
    end
end
