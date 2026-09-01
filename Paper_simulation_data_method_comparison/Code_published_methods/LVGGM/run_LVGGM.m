%% LV-GGM parameter selection and final model fitting
% The script selects lambda and gamma by validation Gaussian negative
% log-likelihood, refits the selected model using all cells, and saves only
% the sparse observed precision matrix S and the low-rank latent-effect
% matrix L from the decomposition K = S - L.

clear;
clc;

%% Analysis settings

% Paths are resolved relative to the directory containing this .m file.
input_relative_path = 'Data_whole_rep1.txt';
logdetppa_relative_path = 'LogdetPPA-0';
output_folder_name = 'LVGGM_Output_rep1';

% Expected number of observed local features in the input matrix.
expected_n_genes = 500;

% Train/validation split used for parameter selection.
rng_seed = 1;
train_fraction = 0.80;

% Candidate regularization parameters.
lambda_grid = [0.02, 0.05, 0.10, 0.20];
gamma_grid  = [2, 5, 10, 20];

% LogdetPPA block parameters.
mu_K = 1;
mu_aux = 1e-8;

% LogdetPPA numerical settings.
OPTIONS.sigma = 10;
OPTIONS.scale_data = 2;
OPTIONS.plotyes = 0;
OPTIONS.tol = 1e-6;
OPTIONS.smoothing = 1;
OPTIONS.printlevel = 0;
OPTIONS.maxiter = 100;

%% Resolve project paths

script_file = mfilename('fullpath');
if isempty(script_file)
    script_dir = pwd;
else
    script_dir = fileparts(script_file);
end

input_file = fullfile(script_dir, input_relative_path);
LogdetPPA_dir = fullfile(script_dir, logdetppa_relative_path);
output_dir = fullfile(script_dir, output_folder_name);

if ~isfolder(LogdetPPA_dir)
    error('LogdetPPA directory not found: %s', LogdetPPA_dir);
end

if ~isfile(input_file)
    error('Input data file not found: %s', input_file);
end

% Initialize a clean output directory for the final LV-GGM matrices.
if isfolder(output_dir)
    rmdir(output_dir, 's');
end
mkdir(output_dir);

addpath(genpath(LogdetPPA_dir));
rehash toolboxcache;

if exist('logdetPPA', 'file') ~= 2
    error('Cannot find logdetPPA.m under: %s', LogdetPPA_dir);
end

%% Read and validate the observed data matrix

X = readmatrix(input_file);

if isempty(X) || ~isnumeric(X)
    error('Input file must contain a numeric matrix.');
end

if any(~isfinite(X(:)))
    error('Input data contain NaN or Inf.');
end

[p, n_cells] = size(X);

if p ~= expected_n_genes
    error('Expected %d genes, but the input contains %d rows.', ...
        expected_n_genes, p);
end

row_sd = std(X, 0, 2);
if any(row_sd < 1e-12)
    bad_rows = find(row_sd < 1e-12, 1, 'first');
    error('Input contains a nearly constant gene row: %d', bad_rows);
end

% Remove only residual row-mean offsets. The input scaling is otherwise
% preserved.
X = X - mean(X, 2);

%% Construct train, validation, and complete covariance matrices

if train_fraction <= 0 || train_fraction >= 1
    error('train_fraction must be strictly between 0 and 1.');
end

rng(rng_seed);
perm = randperm(n_cells);
n_train = floor(train_fraction * n_cells);

if n_train < 1 || n_train >= n_cells
    error('The train/validation split does not contain both subsets.');
end

train_idx = perm(1:n_train);
val_idx = perm(n_train + 1:end);

X_train = X(:, train_idx);
X_val = X(:, val_idx);

Sigma_train = (X_train * X_train') / size(X_train, 2);
Sigma_train = 0.5 * (Sigma_train + Sigma_train');

Sigma_val = (X_val * X_val') / size(X_val, 2);
Sigma_val = 0.5 * (Sigma_val + Sigma_val');

Sigma_all = (X * X') / n_cells;
Sigma_all = 0.5 * (Sigma_all + Sigma_all');

%% Select lambda and gamma by validation negative log-likelihood

validation_nll = inf(numel(lambda_grid), numel(gamma_grid));

for i = 1:numel(lambda_grid)
    for j = 1:numel(gamma_grid)
        lambda = lambda_grid(i);
        gamma = gamma_grid(j);

        try
            [~, ~, K_hat, info] = solve_lvggm_logdetppa( ...
                Sigma_train, lambda, gamma, mu_K, mu_aux, OPTIONS);

            [termcode, pinfeas, dinfeas] = parse_logdetppa_info(info);
            candidate_nll = gaussian_nll_from_precision(Sigma_val, K_hat);

            solver_valid = ...
                (isnan(termcode) || termcode == 0 || termcode == -1) && ...
                (isnan(pinfeas) || pinfeas <= 1e-5) && ...
                (isnan(dinfeas) || dinfeas <= 1e-5);

            if solver_valid && isfinite(candidate_nll)
                validation_nll(i, j) = candidate_nll;
            end
        catch
            validation_nll(i, j) = Inf;
        end
    end
end

[best_nll, best_linear_index] = min(validation_nll(:));

if ~isfinite(best_nll)
    error('All candidate LV-GGM models failed or produced invalid validation likelihoods.');
end

[best_lambda_index, best_gamma_index] = ind2sub( ...
    size(validation_nll), best_linear_index);

lambda_selected = lambda_grid(best_lambda_index);
gamma_selected = gamma_grid(best_gamma_index);

%% Refit the selected model using all cells

[S_final, L_final] = solve_lvggm_logdetppa( ...
    Sigma_all, lambda_selected, gamma_selected, mu_K, mu_aux, OPTIONS);

S_final = 0.5 * (S_final + S_final');
L_final = 0.5 * (L_final + L_final');

%% Save final LV-GGM matrices

% LVGGM_S_local_local_precision.txt
% Complete sparse conditional precision matrix S among observed local
% features. Diagonal entries and unthresholded off-diagonal coefficients are
% retained.
write_matrix_high_precision(fullfile(output_dir, ...
    'LVGGM_S_local_local_precision.txt'), S_final);

% LVGGM_L_lowrank_latent_effect_matrix.txt
% Symmetric low-rank positive-semidefinite latent-effect matrix L in the
% decomposition K = S - L.
write_matrix_high_precision(fullfile(output_dir, ...
    'LVGGM_L_lowrank_latent_effect_matrix.txt'), L_final);

%% Local functions

function [S_hat, L_hat, K_hat, info] = solve_lvggm_logdetppa( ...
    Sigma, lambda, gamma, mu_K, mu_aux, OPTIONS)

    p = size(Sigma, 1);
    n2 = p * (p + 1) / 2;

    % Full-entry L1 weights under the symmetric-vectorization convention.
    l1_weight = zeros(n2, 1);
    pos = 0;

    for col = 1:p
        for row = 1:col
            pos = pos + 1;
            if row == col
                l1_weight(pos) = 1;
            else
                l1_weight(pos) = sqrt(2);
            end
        end
    end

    blk = cell(3, 2);
    At = cell(3, 1);
    C = cell(3, 1);

    % K is the marginal precision matrix and forms the log-determinant block.
    blk{1, 1} = 's';
    blk{1, 2} = p;
    C{1, 1} = Sigma;

    % L is the positive-semidefinite low-rank latent-effect matrix.
    blk{2, 1} = 's';
    blk{2, 2} = p;
    C{2, 1} = lambda * eye(p);

    % The positive vectors U and V represent the absolute-value penalty on S.
    blk{3, 1} = 'l';
    blk{3, 2} = 2 * n2;
    C{3, 1} = lambda * gamma * [l1_weight; l1_weight];

    % The equality constraint implements S = K + L and svec(S) = U - V.
    I_n2 = speye(n2);
    At{1, 1} = I_n2;
    At{2, 1} = I_n2;
    At{3, 1} = [-I_n2, I_n2]';

    b = zeros(n2, 1);
    mu = [mu_K; mu_aux; mu_aux];

    [~, Xsol, ~, ~, info, ~] = logdetPPA( ...
        blk, At, C, b, mu, OPTIONS);

    K_hat = full(Xsol{1});
    K_hat = 0.5 * (K_hat + K_hat');

    L_hat = full(Xsol{2});
    L_hat = 0.5 * (L_hat + L_hat');

    S_hat = K_hat + L_hat;
    S_hat = 0.5 * (S_hat + S_hat');
end

function nll = gaussian_nll_from_precision(Sigma, K)
    K = 0.5 * (K + K');

    [R, flag] = chol(K);
    if flag ~= 0
        nll = Inf;
        return;
    end

    logdetK = 2 * sum(log(diag(R)));
    nll = trace(Sigma * K) - logdetK;
end

function [termcode, pinfeas, dinfeas] = parse_logdetppa_info(info)
    termcode = NaN;
    pinfeas = NaN;
    dinfeas = NaN;

    if isfield(info, 'termcode')
        termcode = info.termcode;
    end

    if isfield(info, 'pinfeas')
        pinfeas = info.pinfeas;
    end

    if isfield(info, 'dinfeas')
        dinfeas = info.dinfeas;
    end
end

function write_matrix_high_precision(filename, M)
    fid = fopen(filename, 'w');

    if fid == -1
        error('Cannot open output file: %s', filename);
    end

    cleanup_object = onCleanup(@() fclose(fid)); %#ok<NASGU>

    [n_rows, n_cols] = size(M);
    row_format = [repmat('%.16e\t', 1, n_cols - 1), '%.16e\n'];

    for row = 1:n_rows
        fprintf(fid, row_format, M(row, :));
    end
end
