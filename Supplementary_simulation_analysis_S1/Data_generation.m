clear;
clc;

%% Simulation structure
% LL(local-local): observed block
% Generates matched full-model and LL-only Gaussian datasets.
% The full model contains 500 local nodes, one hidden node, and one global
% node. Its local network contains 500 within-group edges and 100
% between-group edges. The hidden node connects to 75 local nodes.
% The LL-only model uses the final 500 x 500 LL block from the matched full
% model as its complete precision matrix.

script_file = mfilename('fullpath');
if isempty(script_file)
    script_dir = pwd;
else
    script_dir = fileparts(script_file);
end

RootOutputFolder = fullfile(script_dir, ...
    'Simulated_full_and_LLonly');
reset_output_folder = true;

Alpha_list_file = fullfile(script_dir, 'Alpha_list.txt');
Lambda_list_file = fullfile(script_dir, 'Lambda_list.txt');

if ~isfile(Alpha_list_file)
    error('Cannot find Alpha_list.txt beside the MATLAB script.');
end
if ~isfile(Lambda_list_file)
    error('Cannot find Lambda_list.txt beside the MATLAB script.');
end

if reset_output_folder && isfolder(RootOutputFolder)
    rmdir(RootOutputFolder, 's');
end
if ~isfolder(RootOutputFolder)
    mkdir(RootOutputFolder);
end

%% Network dimensions
n_nodes_per_group = 5;
n_groups = 100;
n_genes = n_nodes_per_group * n_groups;

n_hidden = 1;
n_global = 1;
n_cells = 50000;
n_repeats = 3;

p_full = n_genes + n_hidden + n_global;
idx_genes = 1:n_genes;
idx_hidden = n_genes + 1;
idx_global = n_genes + 2;

%% Network edges
n_within_LL_edges = n_groups * n_nodes_per_group;
n_between_LL_edges = 100;
n_LH_edges = 75;
n_hidden_target_groups = n_LH_edges / n_nodes_per_group;
n_hidden_idx_per_group = 2;

if mod(n_LH_edges, n_nodes_per_group) ~= 0
    error('n_LH_edges must be divisible by n_nodes_per_group.');
end

noise_variance_fraction_grid = [0, 0.10, 0.25, 0.50];
if any(noise_variance_fraction_grid < 0) || ...
        any(noise_variance_fraction_grid > 1)
    error('Noise variance fractions must lie between zero and one.');
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
numerical_zero_threshold = 1e-10;

%% Local groups and between-group pairs
group_id = repelem((1:n_groups)', n_nodes_per_group);
all_local_pairs = nchoosek(1:n_genes, 2);
is_between_pair = ...
    group_id(all_local_pairs(:, 1)) ~= group_id(all_local_pairs(:, 2));
between_gene_pairs = all_local_pairs(is_between_pair, :);
n_possible_between_edges = size(between_gene_pairs, 1);

if n_between_LL_edges > n_possible_between_edges
    error('The requested between-group edge count is not feasible.');
end

%% Repeat-specific LL and LH structures
between_edge_order_by_rep = cell(n_repeats, 1);
between_edge_weight_by_rep = cell(n_repeats, 1);
LH_group_order_by_rep = cell(n_repeats, 1);
LH_group_weight_by_rep = cell(n_repeats, 1);

for rep = 1:n_repeats
    rng(900000 + rep, 'twister');
    between_edge_order_by_rep{rep} = ...
        randperm(n_possible_between_edges);
    between_edge_weight_by_rep{rep} = ...
        between_sign * (between_w_min + ...
        (between_w_max - between_w_min) * ...
        rand(n_between_LL_edges, 1));

    rng(910000 + rep, 'twister');
    LH_group_order_by_rep{rep} = randperm(n_groups);
    LH_group_weight_by_rep{rep} = ...
        -(hidden_w_min + ...
        (hidden_w_max - hidden_w_min) * ...
        rand(n_groups, n_nodes_per_group));
end

%% Combined summary allocation
n_summary_rows = ...
    2 * n_repeats * numel(noise_variance_fraction_grid);

Condition = strings(n_summary_rows, 1);
OutputFolder = strings(n_summary_rows, 1);
RepeatID = zeros(n_summary_rows, 1);
Requested_between_LL_edges = zeros(n_summary_rows, 1);
Requested_LH_edges = zeros(n_summary_rows, 1);
Selected_hidden_groups = zeros(n_summary_rows, 1);
Idx_h_per_hidden_group = zeros(n_summary_rows, 1);
Actual_LL_within_edges = zeros(n_summary_rows, 1);
Actual_LL_between_edges = zeros(n_summary_rows, 1);
Actual_LL_total_edges = zeros(n_summary_rows, 1);
Actual_LH_edges = zeros(n_summary_rows, 1);
Actual_LG_edges = zeros(n_summary_rows, 1);
Idx_h_count_cpp_0based = zeros(n_summary_rows, 1);
Idx_g_count_cpp_0based = zeros(n_summary_rows, 1);
summary_n_cells = zeros(n_summary_rows, 1);
n_local_genes = zeros(n_summary_rows, 1);
n_total_nodes = zeros(n_summary_rows, 1);
summary_row = 0;

%% Generate matched repeats
for rep = 1:n_repeats

    rng(110100 + rep, 'twister');

    %% Full-model LL block
    Theta_LL_base = zeros(n_genes, n_genes);

    for g = 1:n_groups
        group_idx = ...
            ((g - 1) * n_nodes_per_group + 1):...
            (g * n_nodes_per_group);

        Unit_LL = local_unit_diag * eye(n_nodes_per_group);

        for j = 1:(n_nodes_per_group - 1)
            edge_weight = local_unit_sign * ...
                (local_unit_w_min + ...
                (local_unit_w_max - local_unit_w_min) * rand());
            Unit_LL(j, j + 1) = edge_weight;
            Unit_LL(j + 1, j) = edge_weight;
        end

        edge_weight = local_unit_sign * ...
            (local_unit_w_min + ...
            (local_unit_w_max - local_unit_w_min) * rand());
        Unit_LL(1, n_nodes_per_group) = edge_weight;
        Unit_LL(n_nodes_per_group, 1) = edge_weight;

        Theta_LL_base(group_idx, group_idx) = Unit_LL;
    end

    selected_between_ids = ...
        between_edge_order_by_rep{rep}(1:n_between_LL_edges);
    selected_between_pairs = ...
        between_gene_pairs(selected_between_ids, :);
    selected_between_weights = between_edge_weight_by_rep{rep};

    for e = 1:n_between_LL_edges
        i = selected_between_pairs(e, 1);
        j = selected_between_pairs(e, 2);
        w = selected_between_weights(e);
        Theta_LL_base(i, j) = w;
        Theta_LL_base(j, i) = w;
    end

    %% Full-model LH block
    hidden_target_groups = sort(...
        LH_group_order_by_rep{rep}(1:n_hidden_target_groups));
    Theta_LH_base = zeros(n_genes, 1);
    hidden_targets = zeros(n_LH_edges, 1);
    target_counter = 0;

    for k = 1:n_hidden_target_groups
        g = hidden_target_groups(k);
        group_idx = ...
            ((g - 1) * n_nodes_per_group + 1):...
            (g * n_nodes_per_group);

        Theta_LH_base(group_idx) = ...
            LH_group_weight_by_rep{rep}(g, :)';

        hidden_targets(target_counter + (1:n_nodes_per_group)) = ...
            group_idx(:);
        target_counter = target_counter + n_nodes_per_group;
    end
    hidden_targets = sort(hidden_targets);

    %% Full-model LG block
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

    Theta_LG_base = zeros(n_genes, 1);
    for g = 1:n_groups
        group_idx = ...
            ((g - 1) * n_nodes_per_group + 1):...
            (g * n_nodes_per_group);
        Theta_LG_base(group_idx) = global_repeat_weights';
    end

    %% Full precision and standardized covariance
    Theta_unshifted = zeros(p_full, p_full);
    Theta_unshifted(idx_genes, idx_genes) = Theta_LL_base;
    Theta_unshifted(idx_genes, idx_hidden) = Theta_LH_base;
    Theta_unshifted(idx_hidden, idx_genes) = Theta_LH_base';
    Theta_unshifted(idx_genes, idx_global) = Theta_LG_base;
    Theta_unshifted(idx_global, idx_genes) = Theta_LG_base';
    Theta_unshifted(idx_hidden, idx_hidden) = hidden_diag;
    Theta_unshifted(idx_global, idx_global) = global_diag;
    Theta_unshifted = (Theta_unshifted + Theta_unshifted') / 2;

    minimum_eigenvalue = min(eig(Theta_unshifted));
    diagonal_shift = max(0, spd_margin - minimum_eigenvalue);
    Theta_generated = ...
        Theta_unshifted + diagonal_shift * eye(p_full);
    Theta_generated = (Theta_generated + Theta_generated') / 2;

    [~, generated_flag] = chol(Theta_generated);
    if generated_flag ~= 0
        error('The generated full precision matrix is not SPD.');
    end

    Sigma_raw = Theta_generated \ eye(p_full);
    Sigma_raw = (Sigma_raw + Sigma_raw') / 2;
    sd_vector = sqrt(diag(Sigma_raw));

    if any(~isfinite(sd_vector)) || any(sd_vector <= 0)
        error('The full covariance matrix has an invalid diagonal.');
    end

    Sigma_full = Sigma_raw ./ (sd_vector * sd_vector');
    Sigma_full = (Sigma_full + Sigma_full') / 2;
    Theta_full = Sigma_full \ eye(p_full);
    Theta_full = (Theta_full + Theta_full') / 2;
    Theta_full(abs(Theta_full) < numerical_zero_threshold) = 0;

    [~, full_theta_flag] = chol(Theta_full);
    if full_theta_flag ~= 0
        error('The standardized full precision matrix is not SPD.');
    end

    Theta_LL_final = Theta_full(idx_genes, idx_genes);
    Theta_LL_final = (Theta_LL_final + Theta_LL_final') / 2;

    %% Full-model signal and matched noise
    [L_full, sigma_full_flag] = chol(Sigma_full, 'lower');
    if sigma_full_flag ~= 0
        error('The full covariance matrix is not SPD.');
    end

    Z_clean_full = randn(n_cells, p_full) * L_full';
    Z_signal_full = standardize_columns(Z_clean_full);
    clear Z_clean_full
    Noise_full = standardize_columns(randn(size(Z_signal_full)));

    %% Full-model priors
    Idx_h_list_1based = zeros(...
        n_hidden_target_groups * n_hidden_idx_per_group, 1);
    prior_counter = 0;

    for k = 1:n_hidden_target_groups
        g = hidden_target_groups(k);
        group_idx = ...
            ((g - 1) * n_nodes_per_group + 1):...
            (g * n_nodes_per_group);

        selected_positions = ...
            randperm(n_nodes_per_group, n_hidden_idx_per_group);
        selected_indices = group_idx(selected_positions);
        Idx_h_list_1based(...
            prior_counter + (1:n_hidden_idx_per_group)) = ...
            selected_indices(:);
        prior_counter = prior_counter + n_hidden_idx_per_group;
    end

    Idx_h_list = sort(Idx_h_list_1based - 1);
    Idx_g_list = sort(...
        setdiff((1:n_genes)', Idx_h_list_1based) - 1);

    if numel(Idx_h_list) ~= ...
            n_hidden_target_groups * n_hidden_idx_per_group
        error('The full-model hidden-prior count is incorrect.');
    end
    if numel(Idx_h_list) + numel(Idx_g_list) ~= n_genes
        error('The full-model prior lists do not cover all local nodes.');
    end

    %% Network counts shared by both datasets
    LL_mask = abs(Theta_LL_final) > numerical_zero_threshold;
    LL_mask(1:n_genes+1:end) = false;
    same_group_matrix = group_id == group_id';
    LL_within_mask = LL_mask & same_group_matrix;
    LL_between_mask = LL_mask & ~same_group_matrix;

    actual_LL_within = nnz(triu(LL_within_mask, 1));
    actual_LL_between = nnz(triu(LL_between_mask, 1));
    actual_LL_total = nnz(triu(LL_mask, 1));
    actual_LH = nnz(...
        abs(Theta_full(idx_genes, idx_hidden)) > ...
        numerical_zero_threshold);
    actual_LG = nnz(...
        abs(Theta_full(idx_genes, idx_global)) > ...
        numerical_zero_threshold);

    %% Save full-model datasets
    for iNoise = 1:numel(noise_variance_fraction_grid)
        noise_fraction = noise_variance_fraction_grid(iNoise);
        noise_percent = round(100 * noise_fraction);
        signal_weight = sqrt(1 - noise_fraction);
        noise_weight = sqrt(noise_fraction);

        Z_observed_full = ...
            signal_weight * Z_signal_full + ...
            noise_weight * Noise_full;
        Z_observed_full = standardize_columns(Z_observed_full);
        Matrix_full = Z_observed_full(:, idx_genes)';

        full_condition = sprintf(...
            'Data_LLbetween_%04d_LH_%03d_noise_%03dpct_rep_%02d', ...
            n_between_LL_edges, n_LH_edges, noise_percent, rep);
        full_data_folder = fullfile(...
            RootOutputFolder, full_condition, 'Data');
        if ~isfolder(full_data_folder)
            mkdir(full_data_folder);
        end

        save_six_outputs(...
            full_data_folder, Theta_full, Matrix_full, ...
            Idx_h_list, Idx_g_list, ...
            Alpha_list_file, Lambda_list_file);

        summary_row = summary_row + 1;
        Condition(summary_row) = string(full_condition);
        OutputFolder(summary_row) = string(full_condition);
        RepeatID(summary_row) = rep;
        Requested_between_LL_edges(summary_row) = n_between_LL_edges;
        Requested_LH_edges(summary_row) = n_LH_edges;
        Selected_hidden_groups(summary_row) = n_hidden_target_groups;
        Idx_h_per_hidden_group(summary_row) = n_hidden_idx_per_group;
        Actual_LL_within_edges(summary_row) = actual_LL_within;
        Actual_LL_between_edges(summary_row) = actual_LL_between;
        Actual_LL_total_edges(summary_row) = actual_LL_total;
        Actual_LH_edges(summary_row) = actual_LH;
        Actual_LG_edges(summary_row) = actual_LG;
        Idx_h_count_cpp_0based(summary_row) = numel(Idx_h_list);
        Idx_g_count_cpp_0based(summary_row) = numel(Idx_g_list);
        summary_n_cells(summary_row) = n_cells;
        n_local_genes(summary_row) = n_genes;
        n_total_nodes(summary_row) = p_full;

        clear Z_observed_full Matrix_full
    end

    clear Z_signal_full Noise_full L_full

    %% LL-only precision and matched data
    [~, ll_theta_flag] = chol(Theta_LL_final);
    if ll_theta_flag ~= 0
        error('The LL block extracted from the full model is not SPD.');
    end

    Sigma_LLonly = Theta_LL_final \ eye(n_genes);
    Sigma_LLonly = (Sigma_LLonly + Sigma_LLonly') / 2;

    [L_LLonly, sigma_ll_flag] = chol(Sigma_LLonly, 'lower');
    if sigma_ll_flag ~= 0
        error('The LL-only covariance matrix is not SPD.');
    end

    rng(200000 + rep, 'twister');
    Z_clean_LLonly = randn(n_cells, n_genes) * L_LLonly';

    rng(300000 + rep, 'twister');
    Noise_LLonly = standardize_columns(randn(n_cells, n_genes));

    Idx_h_LLonly = [];
    Idx_g_LLonly = (0:(n_genes - 1))';

    %% Save LL-only datasets
    for iNoise = 1:numel(noise_variance_fraction_grid)
        noise_fraction = noise_variance_fraction_grid(iNoise);
        noise_percent = round(100 * noise_fraction);
        signal_weight = sqrt(1 - noise_fraction);
        noise_weight = sqrt(noise_fraction);

        Z_observed_LLonly = ...
            signal_weight * Z_clean_LLonly + ...
            noise_weight * Noise_LLonly;
        Matrix_LLonly = Z_observed_LLonly';

        llonly_condition = sprintf(...
            'Data_LLonly_fromFullLL_noise_%03dpct_rep_%02d', ...
            noise_percent, rep);
        llonly_data_folder = fullfile(...
            RootOutputFolder, llonly_condition, 'Data');
        if ~isfolder(llonly_data_folder)
            mkdir(llonly_data_folder);
        end

        save_six_outputs(...
            llonly_data_folder, Theta_LL_final, Matrix_LLonly, ...
            Idx_h_LLonly, Idx_g_LLonly, ...
            Alpha_list_file, Lambda_list_file);

        summary_row = summary_row + 1;
        Condition(summary_row) = string(llonly_condition);
        OutputFolder(summary_row) = string(llonly_condition);
        RepeatID(summary_row) = rep;
        Requested_between_LL_edges(summary_row) = n_between_LL_edges;
        Requested_LH_edges(summary_row) = 0;
        Selected_hidden_groups(summary_row) = 0;
        Idx_h_per_hidden_group(summary_row) = 0;
        Actual_LL_within_edges(summary_row) = actual_LL_within;
        Actual_LL_between_edges(summary_row) = actual_LL_between;
        Actual_LL_total_edges(summary_row) = actual_LL_total;
        Actual_LH_edges(summary_row) = 0;
        Actual_LG_edges(summary_row) = 0;
        Idx_h_count_cpp_0based(summary_row) = 0;
        Idx_g_count_cpp_0based(summary_row) = n_genes;
        summary_n_cells(summary_row) = n_cells;
        n_local_genes(summary_row) = n_genes;
        n_total_nodes(summary_row) = n_genes;

        clear Z_observed_LLonly Matrix_LLonly
    end

    clear Z_clean_LLonly Noise_LLonly Sigma_LLonly L_LLonly
end


%% Save the combined condition summary
All_conditions_summary_Simulated = table(...
    Condition, ...
    OutputFolder, ...
    RepeatID, ...
    Requested_between_LL_edges, ...
    Requested_LH_edges, ...
    Selected_hidden_groups, ...
    Idx_h_per_hidden_group, ...
    Actual_LL_within_edges, ...
    Actual_LL_between_edges, ...
    Actual_LL_total_edges, ...
    Actual_LH_edges, ...
    Actual_LG_edges, ...
    Idx_h_count_cpp_0based, ...
    Idx_g_count_cpp_0based, ...
    summary_n_cells, ...
    n_local_genes, ...
    n_total_nodes, ...
    'VariableNames', {...
        'Condition', ...
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
        'n_cells', ...
        'n_local_genes', ...
        'n_total_nodes'});

writetable(All_conditions_summary_Simulated, ...
    fullfile(RootOutputFolder, ...
    'All_conditions_summary_Simulated.csv'));

%% Local functions
function Xz = standardize_columns(X)
% Centers and scales each column using its sample standard deviation.
    column_mean = mean(X, 1);
    column_sd = std(X, 0, 1);
    column_sd(~isfinite(column_sd) | column_sd == 0) = 1;
    Xz = (X - column_mean) ./ column_sd;

    if any(~isfinite(Xz), 'all')
        error('A standardized matrix contains NaN or Inf.');
    end
end

function save_six_outputs(...
    data_folder, Theta_full, Matrix, Idx_h_list, Idx_g_list, ...
    alpha_file, lambda_file)
% Saves the six files required by the downstream model-fitting workflow.
    writematrix(Theta_full, ...
        fullfile(data_folder, 'Theta_full.txt'), ...
        'Delimiter', 'tab');

    write_matrix_or_empty(...
        fullfile(data_folder, 'Idx_h_list.txt'), Idx_h_list);
    write_matrix_or_empty(...
        fullfile(data_folder, 'Idx_g_list.txt'), Idx_g_list);

    write_binary_matrix(...
        fullfile(data_folder, 'Data_whole.bin'), Matrix);

    [status_alpha, message_alpha] = copyfile(...
        alpha_file, fullfile(data_folder, 'Alpha_list.txt'));
    if ~status_alpha
        error('Unable to copy Alpha_list.txt: %s', message_alpha);
    end

    [status_lambda, message_lambda] = copyfile(...
        lambda_file, fullfile(data_folder, 'Lambda_list.txt'));
    if ~status_lambda
        error('Unable to copy Lambda_list.txt: %s', message_lambda);
    end
end

function write_binary_matrix(bin_file, Matrix)
% Writes a genes-by-samples matrix in the C++ row-major binary layout.
    fid = fopen(bin_file, 'wb');
    if fid == -1
        error('Cannot open binary output file: %s', bin_file);
    end

    cleanup_object = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fwrite(fid, uint64(size(Matrix, 1)), 'uint64');
    fwrite(fid, uint64(size(Matrix, 2)), 'uint64');
    fwrite(fid, Matrix', 'double');
end

function write_matrix_or_empty(filename, X)
% Creates an empty file for an empty index list.
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
