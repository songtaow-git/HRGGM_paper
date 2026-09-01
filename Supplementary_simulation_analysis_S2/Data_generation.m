clear;
clc;
close all;

%% Simulation structure
% Generates Gaussian datasets from precision matrices containing 200 local
% nodes （observed nodes）, one hidden node, and one global node.
%
% The local(observed) network contains 40 groups of five nodes. Each group contributes
% a five-edge ring. The local-hidden block is fixed at 50 edges, while the
% number of between-group local-local edges varies across conditions.
%
% Each condition is generated for three repeats and three observation-noise
% variance fractions. Within a condition and repeat, the precision matrix,
% clean samples, and noise realization are shared across noise levels.
%
% Each condition folder contains a Data subfolder with:
%   Alpha_list.txt
%   Data_whole.bin
%   Idx_g_list.txt
%   Idx_h_list.txt
%   Lambda_list.txt
%   Theta_full.txt

%% Paths
script_file = mfilename('fullpath');

if isempty(script_file)
    base_dir = pwd;
else
    base_dir = fileparts(script_file);
end

RootOutputFolder = fullfile( ...
    base_dir, ...
    'LLdensity_experiment');

Alpha_list_file = fullfile(base_dir, 'Alpha_list.txt');
Lambda_list_file = fullfile(base_dir, 'Lambda_list.txt');

if ~isfile(Alpha_list_file)
    error('Required file not found: %s', Alpha_list_file);
end

if ~isfile(Lambda_list_file)
    error('Required file not found: %s', Lambda_list_file);
end

reset_output_folder = true;

if reset_output_folder && isfolder(RootOutputFolder)
    rmdir(RootOutputFolder, 's');
end

if ~isfolder(RootOutputFolder)
    mkdir(RootOutputFolder);
end

%% Node dimensions
n_nodes_per_group = 5;
n_groups = 40;
n_genes = n_nodes_per_group * n_groups;

n_cells = 10000;
n_hidden = 1;
n_global = 1;
n_repeats = 3;

if n_hidden ~= 1 || n_global ~= 1
    error('The simulation requires one hidden node and one global node.');
end

p = n_genes + n_hidden + n_global;

idx_genes = 1:n_genes;
idx_hidden = n_genes + 1;
idx_global = n_genes + 2;

%% Local-local density conditions
condition_table = table( ...
    ["LL_0000"; "LL_0050"; "LL_0200"; "LL_1000"], ...
    [0; 50; 200; 1000], ...
    repmat("LL_density", 4, 1), ...
    'VariableNames', { ...
        'ConditionLabel', ...
        'Between_LL_edges', ...
        'ExperimentType'});

fixed_n_LH_edges = 50;
n_hidden_target_groups = fixed_n_LH_edges / n_nodes_per_group;

if mod(fixed_n_LH_edges, n_nodes_per_group) ~= 0 || ...
        fixed_n_LH_edges > n_genes
    error('The fixed LH edge count is incompatible with the node structure.');
end

n_conditions = height(condition_table);
noise_variance_fraction_grid = [0, 0.10, 0.25];

if any(noise_variance_fraction_grid < 0) || ...
        any(noise_variance_fraction_grid > 1)
    error('Noise variance fractions must be between zero and one.');
end

%% Precision-matrix parameters
local_unit_diag = 2;
local_unit_w_min = 0.4;
local_unit_w_max = 0.8;
local_unit_sign = -1;

between_w_min = local_unit_w_min;
between_w_max = local_unit_w_max;
between_sign = local_unit_sign;

hidden_diag = 4;
hidden_w_min = 0.45;
hidden_w_max = 0.75;

global_diag = 8;
global_shift = 0.1;
global_span = 0.5;

spd_margin = 0.1;
standardize_cov_diag = true;
numerical_zero_threshold = 1e-10;
add_noise = true;

n_hidden_idx_per_group = 2;

if n_hidden_idx_per_group < 0 || ...
        n_hidden_idx_per_group > n_nodes_per_group
    error('Invalid number of hidden-prior indices per group.');
end

%% Local-node groups and between-group pairs
group_id = repelem((1:n_groups)', n_nodes_per_group);

all_local_gene_pairs = nchoosek(1:n_genes, 2);
is_between_gene_pair = ...
    group_id(all_local_gene_pairs(:, 1)) ~= ...
    group_id(all_local_gene_pairs(:, 2));

between_gene_pairs = all_local_gene_pairs(is_between_gene_pair, :);
n_possible_between_LL_edges = size(between_gene_pairs, 1);

max_requested_between_edges = max(condition_table.Between_LL_edges);

if max_requested_between_edges > n_possible_between_LL_edges
    error('Requested between-group LL edges exceed the available pairs.');
end

%% Repeat-specific nested structures
between_edge_order_by_rep = cell(n_repeats, 1);
between_edge_weight_by_rep = cell(n_repeats, 1);
LH_group_order_by_rep = cell(n_repeats, 1);
LH_group_weight_by_rep = cell(n_repeats, 1);

for rep = 1:n_repeats
    rng(900000 + rep, 'twister');

    between_edge_order_by_rep{rep} = ...
        randperm(n_possible_between_LL_edges);

    between_edge_weight_by_rep{rep} = ...
        between_sign * ...
        (between_w_min + ...
        (between_w_max - between_w_min) * ...
        rand(max_requested_between_edges, 1));

    rng(910000 + rep, 'twister');

    LH_group_order_by_rep{rep} = randperm(n_groups);

    LH_group_weight_by_rep{rep} = ...
        -(hidden_w_min + ...
        (hidden_w_max - hidden_w_min) * ...
        rand(n_groups, n_nodes_per_group));
end

overall_summary_initialized = false;

%% Generate local-local density conditions
for iCond = 1:n_conditions
    condition_label = char(condition_table.ConditionLabel(iCond));
    experiment_type = char(condition_table.ExperimentType(iCond));
    n_requested_between_edges = condition_table.Between_LL_edges(iCond);
    n_LH_edges = fixed_n_LH_edges;

    if n_requested_between_edges < 0 || ...
            n_requested_between_edges > n_possible_between_LL_edges
        error('Invalid between-group LL edge count.');
    end

    if n_hidden_target_groups > n_groups
        error('The number of hidden-connected groups exceeds n_groups.');
    end

    for rep = 1:n_repeats
        rng(100000 + iCond * 1000 + rep, 'twister');

        %% Local-local precision block
        Theta_LL = zeros(n_genes, n_genes);

        for g = 1:n_groups
            group_idx = ...
                ((g - 1) * n_nodes_per_group + 1): ...
                (g * n_nodes_per_group);

            Unit_precision_LL = ...
                local_unit_diag * eye(n_nodes_per_group);

            for j = 1:(n_nodes_per_group - 1)
                edge_weight = ...
                    local_unit_sign * ...
                    (local_unit_w_min + ...
                    (local_unit_w_max - local_unit_w_min) * rand());

                Unit_precision_LL(j, j + 1) = edge_weight;
                Unit_precision_LL(j + 1, j) = edge_weight;
            end

            edge_weight = ...
                local_unit_sign * ...
                (local_unit_w_min + ...
                (local_unit_w_max - local_unit_w_min) * rand());

            Unit_precision_LL(1, n_nodes_per_group) = edge_weight;
            Unit_precision_LL(n_nodes_per_group, 1) = edge_weight;

            Theta_LL(group_idx, group_idx) = Unit_precision_LL;
        end

        if n_requested_between_edges > 0
            selected_between_ids = ...
                between_edge_order_by_rep{rep}( ...
                1:n_requested_between_edges);

            selected_between_pairs = ...
                between_gene_pairs(selected_between_ids, :);

            selected_between_weights = ...
                between_edge_weight_by_rep{rep}( ...
                1:n_requested_between_edges);

            for e = 1:n_requested_between_edges
                gene_i = selected_between_pairs(e, 1);
                gene_j = selected_between_pairs(e, 2);
                edge_weight = selected_between_weights(e);

                Theta_LL(gene_i, gene_j) = edge_weight;
                Theta_LL(gene_j, gene_i) = edge_weight;
            end
        end

        %% Local-hidden precision block
        Theta_LH = zeros(n_genes, 1);

        if n_hidden_target_groups > 0
            hidden_target_groups = sort( ...
                LH_group_order_by_rep{rep}( ...
                1:n_hidden_target_groups));
        else
            hidden_target_groups = [];
        end

        for k = 1:numel(hidden_target_groups)
            g = hidden_target_groups(k);

            group_idx = ...
                ((g - 1) * n_nodes_per_group + 1): ...
                (g * n_nodes_per_group);

            Theta_LH(group_idx, 1) = ...
                LH_group_weight_by_rep{rep}(g, :).';
        end

        %% Local-global precision block
        global_repeat_weights = ...
            (rand(1, n_nodes_per_group) - 0.5) * global_span;

        for j = 1:n_nodes_per_group
            if global_repeat_weights(j) < 0
                global_repeat_weights(j) = ...
                    global_repeat_weights(j) - global_shift;
            elseif global_repeat_weights(j) > 0
                global_repeat_weights(j) = ...
                    global_repeat_weights(j) + global_shift;
            end
        end

        Theta_LG = zeros(n_genes, 1);

        for g = 1:n_groups
            group_idx = ...
                ((g - 1) * n_nodes_per_group + 1): ...
                (g * n_nodes_per_group);

            Theta_LG(group_idx, 1) = global_repeat_weights.';
        end

        %% Complete precision matrix
        Theta_unshifted = zeros(p, p);

        Theta_unshifted(idx_genes, idx_genes) = Theta_LL;

        Theta_unshifted(idx_genes, idx_hidden) = Theta_LH;
        Theta_unshifted(idx_hidden, idx_genes) = Theta_LH.';

        Theta_unshifted(idx_genes, idx_global) = Theta_LG;
        Theta_unshifted(idx_global, idx_genes) = Theta_LG.';

        Theta_unshifted(idx_hidden, idx_hidden) = hidden_diag;
        Theta_unshifted(idx_global, idx_global) = global_diag;

        Theta_unshifted(idx_hidden, idx_global) = 0;
        Theta_unshifted(idx_global, idx_hidden) = 0;

        Theta_unshifted = ...
            0.5 * (Theta_unshifted + Theta_unshifted.');

        minimum_eigenvalue = min(eig(Theta_unshifted));
        diagonal_shift = max(0, spd_margin - minimum_eigenvalue);

        Theta_generated = ...
            Theta_unshifted + diagonal_shift * eye(p);

        Theta_generated = ...
            0.5 * (Theta_generated + Theta_generated.');

        [~, generated_flag] = chol(Theta_generated);

        if generated_flag ~= 0
            error('The shifted precision matrix is not positive definite.');
        end

        %% Covariance standardization
        Sigma_raw = Theta_generated \ eye(p);
        Sigma_raw = 0.5 * (Sigma_raw + Sigma_raw.');

        if standardize_cov_diag
            sd_vec = sqrt(diag(Sigma_raw));

            if any(~isfinite(sd_vec)) || any(sd_vec <= 0)
                error('The covariance matrix has an invalid diagonal.');
            end

            Sigma = Sigma_raw ./ (sd_vec * sd_vec.');
            Sigma = 0.5 * (Sigma + Sigma.');
        else
            Sigma = Sigma_raw;
        end

        Theta = Sigma \ eye(p);
        Theta = 0.5 * (Theta + Theta.');
        Theta(abs(Theta) < numerical_zero_threshold) = 0;

        [~, final_flag] = chol(Theta);

        if final_flag ~= 0
            error('The standardized precision matrix is not positive definite.');
        end

        %% Clean samples and shared noise
        [L_sigma, sigma_flag] = chol(Sigma, 'lower');

        if sigma_flag ~= 0
            error('The covariance matrix is not positive definite.');
        end

        Z_clean = randn(n_cells, p) * L_sigma.';

        signal_mean = mean(Z_clean, 1);
        signal_sd = std(Z_clean, 0, 1);
        signal_sd(signal_sd == 0) = 1;
        Z_signal = (Z_clean - signal_mean) ./ signal_sd;

        if add_noise
            Noise_raw = randn(size(Z_signal));
            noise_mean = mean(Noise_raw, 1);
            noise_sd = std(Noise_raw, 0, 1);
            noise_sd(noise_sd == 0) = 1;
            Noise_z = (Noise_raw - noise_mean) ./ noise_sd;
        else
            Noise_z = zeros(size(Z_signal));
        end

        %% Prior and remaining local-node indices
        Idx_h_list_1based = [];

        for k = 1:numel(hidden_target_groups)
            g = hidden_target_groups(k);

            group_idx_1based = ...
                ((g - 1) * n_nodes_per_group + 1): ...
                (g * n_nodes_per_group);

            selected_order = randperm( ...
                n_nodes_per_group, ...
                n_hidden_idx_per_group);

            selected_idx_1based = ...
                group_idx_1based(selected_order);

            Idx_h_list_1based = [ ...
                Idx_h_list_1based; ...
                selected_idx_1based(:)]; %#ok<AGROW>
        end

        all_local_idx_1based = (1:n_genes).';

        Idx_g_list_1based = setdiff( ...
            all_local_idx_1based, ...
            Idx_h_list_1based);

        Idx_h_list = sort(Idx_h_list_1based(:) - 1);
        Idx_g_list = sort(Idx_g_list_1based(:) - 1);

        expected_n_h = ...
            numel(hidden_target_groups) * ...
            n_hidden_idx_per_group;

        if numel(Idx_h_list) ~= expected_n_h
            error('Unexpected number of hidden-prior indices.');
        end

        if numel(Idx_h_list) + numel(Idx_g_list) ~= n_genes
            error('Idx_h_list and Idx_g_list do not cover all local nodes.');
        end

        if ~isempty(intersect(Idx_h_list, Idx_g_list))
            error('Idx_h_list and Idx_g_list overlap.');
        end

        if ~isempty(Idx_h_list) && ...
                (min(Idx_h_list) < 0 || max(Idx_h_list) > n_genes - 1)
            error('Idx_h_list contains an invalid C++ index.');
        end

        if ~isempty(Idx_g_list) && ...
                (min(Idx_g_list) < 0 || max(Idx_g_list) > n_genes - 1)
            error('Idx_g_list contains an invalid C++ index.');
        end

        actual_LL_within_edges = ...
            n_groups * n_nodes_per_group;

        actual_LL_between_edges = n_requested_between_edges;

        actual_LL_total_edges = ...
            actual_LL_within_edges + ...
            actual_LL_between_edges;

        actual_LH_edges = nnz(Theta_LH);
        actual_LG_edges = nnz(Theta_LG);

        %% Noise-level datasets
        for iNoise = 1:numel(noise_variance_fraction_grid)
            noise_variance_fraction = ...
                noise_variance_fraction_grid(iNoise);

            noise_percent = ...
                round(noise_variance_fraction * 100);

            Z_observed = ...
                sqrt(1 - noise_variance_fraction) * Z_signal + ...
                sqrt(noise_variance_fraction) * Noise_z;

            observed_mean = mean(Z_observed, 1);
            observed_sd = std(Z_observed, 0, 1);
            observed_sd(observed_sd == 0) = 1;

            Z_observed = ...
                (Z_observed - observed_mean) ./ observed_sd;

            X_local = Z_observed(:, idx_genes).';

            condition_name = sprintf( ...
                ['Data_%s_LLbetween_%04d_LH_%03d_' ...
                 'noise_%03dpct_rep_%02d'], ...
                condition_label, ...
                n_requested_between_edges, ...
                n_LH_edges, ...
                noise_percent, ...
                rep);

            condition_folder = fullfile( ...
                RootOutputFolder, ...
                condition_name);

            data_folder = fullfile(condition_folder, 'Data');

            if ~isfolder(data_folder)
                mkdir(data_folder);
            end

            writematrix( ...
                Theta, ...
                fullfile(data_folder, 'Theta_full.txt'), ...
                'Delimiter', 'tab');

            write_matrix_or_empty( ...
                fullfile(data_folder, 'Idx_h_list.txt'), ...
                Idx_h_list);

            write_matrix_or_empty( ...
                fullfile(data_folder, 'Idx_g_list.txt'), ...
                Idx_g_list);

            copyfile( ...
                Alpha_list_file, ...
                fullfile(data_folder, 'Alpha_list.txt'));

            copyfile( ...
                Lambda_list_file, ...
                fullfile(data_folder, 'Lambda_list.txt'));

            write_binary_matrix( ...
                fullfile(data_folder, 'Data_whole.bin'), ...
                X_local);

            condition_summary = table( ...
                string(condition_name), ...
                string(condition_label), ...
                string(experiment_type), ...
                string(condition_folder), ...
                rep, ...
                n_requested_between_edges, ...
                n_LH_edges, ...
                n_hidden_target_groups, ...
                n_hidden_idx_per_group, ...
                actual_LL_within_edges, ...
                actual_LL_between_edges, ...
                actual_LL_total_edges, ...
                actual_LH_edges, ...
                actual_LG_edges, ...
                numel(Idx_h_list), ...
                numel(Idx_g_list), ...
                noise_variance_fraction, ...
                n_cells, ...
                n_genes, ...
                p, ...
                'VariableNames', { ...
                    'Condition', ...
                    'ConditionLabel', ...
                    'ExperimentType', ...
                    'OutputFolder', ...
                    'RepeatID', ...
                    'Requested_between_LL_edges', ...
                    'Requested_LH_edges', ...
                    'Selected_hidden_groups', ...
                    'Idx_h_per_hidden_group', ...
                    'Actual_LL_within_edges', ...
                    'Actual_LL_between_edges', ...
                    'Actual_LL_total_edges', ...
                    'Actual_LH_edges', ...
                    'Actual_LG_edges', ...
                    'Idx_h_count_cpp_0based', ...
                    'Idx_g_count_cpp_0based', ...
                    'Noise_variance_fraction', ...
                    'n_cells', ...
                    'n_local_genes', ...
                    'n_total_nodes'});

            if ~overall_summary_initialized
                overall_summary = condition_summary;
                overall_summary_initialized = true;
            else
                overall_summary = [ ...
                    overall_summary; ...
                    condition_summary]; %#ok<AGROW>
            end

            clear Z_observed X_local
        end

        clear Z_clean Z_signal Noise_raw Noise_z
        clear Sigma_raw Sigma Theta Theta_generated Theta_unshifted
    end
end

if overall_summary_initialized
    writetable( ...
        overall_summary, ...
        fullfile(RootOutputFolder, 'All_conditions_summary.csv'));
end

%% Local functions
function write_binary_matrix(bin_file, Matrix)
% Writes a matrix using uint64 dimensions followed by row-major doubles.

    fid = fopen(bin_file, 'wb');

    if fid == -1
        error('Cannot open binary output file: %s', bin_file);
    end

    cleanup_object = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fwrite(fid, uint64(size(Matrix, 1)), 'uint64');
    fwrite(fid, uint64(size(Matrix, 2)), 'uint64');
    fwrite(fid, Matrix.', 'double');
end

function write_matrix_or_empty(filename, X)
% Writes a tab-delimited vector or creates an empty file.

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