%% Network recovery among observed variables
% Compares six methods across three matched repeats using the upper
% triangle of the symmetrized 500-by-500 observed-variable precision block.
% Outputs repeat-level and summary metrics, mean precision-recall curves,
% candidate-list recovery results, and a two-panel publication figure.

clear;
clc;
close all;

% Set default fonts for generated figures.
set(groot, ...
    'defaultAxesFontName', 'Arial', ...
    'defaultTextFontName', 'Arial', ...
    'defaultLegendFontName', 'Arial');

%% Input, output, and analysis settings
script_file = mfilename('fullpath');

if isempty(script_file)
    root_dir = pwd;
else
    root_dir = fileparts(script_file);
end

output_dir = fullfile(root_dir, 'Observed_variable_network_analysis');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

Nlocal = 500;

Nrep = 3;

GT_zero_tolerance = 1e-12;

target_recall = 0.95;

candidate_multipliers = [0.5, 1, 2, 3, 5, 10, 20];

recall_grid = (0:0.001:1)';

%% Method definitions
methods = struct([]);

methods(1).display_name = 'HR-GGM';
methods(1).folder       = 'HRGGM';
methods(1).file_prefix  = 'HRGGM';

methods(2).display_name = 'LV-GGM';
methods(2).folder       = 'LVGGM';
methods(2).file_prefix  = 'LVGGM';

methods(3).display_name = 'TGL';
methods(3).folder       = 'TGL';
methods(3).file_prefix  = 'TGL';

methods(4).display_name = 'GRNBoost2';
methods(4).folder       = 'GRNBoost2';
methods(4).file_prefix  = 'GRNBoost2';

methods(5).display_name = 'SCING';
methods(5).folder       = 'SCING';
methods(5).file_prefix  = 'SCING';

methods(6).display_name = 'RegDiffusion';
methods(6).folder       = 'RegDiffusion';
methods(6).file_prefix  = 'RegDiffusion';

Nmethod = numel(methods);

GT_folder = fullfile(root_dir, 'GT');
GT_prefix = 'GT';

%% Result storage
Nresult = Nmethod * Nrep;

Method             = strings(Nresult, 1);
Repeat             = zeros(Nresult, 1);
TrueEdges          = zeros(Nresult, 1);
TotalPossibleEdges = zeros(Nresult, 1);
PositiveFraction   = zeros(Nresult, 1);

AUPRC           = zeros(Nresult, 1);
NormalizedAUPRC = zeros(Nresult, 1);
TopKRecovery    = zeros(Nresult, 1);

Budget95Edges         = zeros(Nresult, 1);
Budget95Factor        = zeros(Nresult, 1);
PrecisionAt95Recall   = zeros(Nresult, 1);
FalsePositiveBurden95 = zeros(Nresult, 1);

PR_precision_grid = nan( ...
    numel(recall_grid), ...
    Nmethod, ...
    Nrep);

candidate_recovery = nan( ...
    numel(candidate_multipliers), ...
    Nmethod, ...
    Nrep);

result_row = 0;

%% Analyze matched repeats
for rep = 1:Nrep

    % Read and preprocess the matched ground-truth matrix.
    GT_file = resolve_matrix_file( ...
        GT_folder, ...
        GT_prefix, ...
        rep);

    GT_full = read_numeric_matrix(GT_file);

    if size(GT_full, 1) < Nlocal || size(GT_full, 2) < Nlocal
        error(['GT repeat %d has size %d x %d, but at least ', ...
               '%d x %d is required.'], ...
               rep, ...
               size(GT_full, 1), ...
               size(GT_full, 2), ...
               Nlocal, ...
               Nlocal);
    end

    % Symmetrize the full matrix before extracting the observed-variable block.
    GT_full = symmetrize_matrix(GT_full);

    GT_LL = GT_full(1:Nlocal, 1:Nlocal);
    GT_LL = abs(GT_LL);

    % Exclude diagonal entries from edge evaluation.
    GT_LL(1:Nlocal+1:end) = 0;

    upper_mask = triu(true(Nlocal), 1);

    truth_labels = GT_LL(upper_mask) > GT_zero_tolerance;
    truth_labels = logical(truth_labels(:));

    K = sum(truth_labels);
    Npairs = numel(truth_labels);
    prevalence = K / Npairs;

    if K == 0
        error('No true LL edges were found in GT repeat %d.', rep);
    end

    % Evaluate each estimated network against the same ground truth.
    for m = 1:Nmethod

        result_row = result_row + 1;

        method_folder = fullfile( ...
            root_dir, ...
            methods(m).folder);

        method_file = resolve_matrix_file( ...
            method_folder, ...
            methods(m).file_prefix, ...
            rep);

        estimated_full = read_numeric_matrix(method_file);

        if size(estimated_full, 1) < Nlocal || ...
                size(estimated_full, 2) < Nlocal

            error(['%s repeat %d has size %d x %d, but at least ', ...
                   '%d x %d is required.'], ...
                   methods(m).display_name, ...
                   rep, ...
                   size(estimated_full, 1), ...
                   size(estimated_full, 2), ...
                   Nlocal, ...
                   Nlocal);
        end

        % Symmetrize the full estimate before extracting the observed-variable block.
        estimated_full = symmetrize_matrix(estimated_full);

        estimated_LL = estimated_full(1:Nlocal, 1:Nlocal);
        estimated_LL = abs(estimated_LL);

        % Exclude diagonal entries from edge evaluation.
        estimated_LL(1:Nlocal+1:end) = 0;

        edge_scores = estimated_LL(upper_mask);
        edge_scores = edge_scores(:);

        % Replace non-finite edge scores with zero before ranking.
        if any(~isfinite(edge_scores))

            edge_scores(~isfinite(edge_scores)) = 0;
        end

        % Calculate tie-aware ranking and candidate-budget metrics.
        rank_result = calculate_ranking_metrics( ...
            edge_scores, ...
            truth_labels, ...
            K, ...
            target_recall, ...
            candidate_multipliers, ...
            recall_grid);

        normalized_auprc = ...
            (rank_result.AUPRC - prevalence) / ...
            max(1 - prevalence, eps);

        Method(result_row)             = methods(m).display_name;
        Repeat(result_row)             = rep;
        TrueEdges(result_row)          = K;
        TotalPossibleEdges(result_row) = Npairs;
        PositiveFraction(result_row)   = prevalence;

        AUPRC(result_row)           = rank_result.AUPRC;
        NormalizedAUPRC(result_row) = normalized_auprc;
        TopKRecovery(result_row)    = rank_result.TopKRecovery;

        Budget95Edges(result_row) = ...
            rank_result.TargetBudgetEdges;

        Budget95Factor(result_row) = ...
            rank_result.TargetBudgetFactor;

        PrecisionAt95Recall(result_row) = ...
            rank_result.PrecisionAtTarget;

        FalsePositiveBurden95(result_row) = ...
            rank_result.FalsePositiveBurden;

        PR_precision_grid(:, m, rep) = ...
            rank_result.PrecisionOnRecallGrid;

        candidate_recovery(:, m, rep) = ...
            rank_result.CandidateRecovery;

    end
end

%% Build numerical output tables
% Store one row for each method and repeat.
per_repeat_table = table( ...
    Method, ...
    Repeat, ...
    TrueEdges, ...
    TotalPossibleEdges, ...
    PositiveFraction, ...
    AUPRC, ...
    NormalizedAUPRC, ...
    TopKRecovery, ...
    Budget95Edges, ...
    Budget95Factor, ...
    PrecisionAt95Recall, ...
    FalsePositiveBurden95);

% Summarize metrics across repeats for each method.
summary_method = strings(Nmethod, 1);

mean_AUPRC = zeros(Nmethod, 1);
sd_AUPRC   = zeros(Nmethod, 1);

mean_nAUPRC = zeros(Nmethod, 1);
sd_nAUPRC   = zeros(Nmethod, 1);

mean_TopK = zeros(Nmethod, 1);
sd_TopK   = zeros(Nmethod, 1);

mean_BudgetFactor = zeros(Nmethod, 1);
sd_BudgetFactor   = zeros(Nmethod, 1);

mean_Precision95 = zeros(Nmethod, 1);
sd_Precision95   = zeros(Nmethod, 1);

for m = 1:Nmethod

    method_name = methods(m).display_name;

    use_rows = per_repeat_table.Method == method_name;

    summary_method(m) = method_name;

    mean_AUPRC(m) = mean( ...
        per_repeat_table.AUPRC(use_rows), ...
        'omitnan');

    sd_AUPRC(m) = std( ...
        per_repeat_table.AUPRC(use_rows), ...
        0, ...
        'omitnan');

    mean_nAUPRC(m) = mean( ...
        per_repeat_table.NormalizedAUPRC(use_rows), ...
        'omitnan');

    sd_nAUPRC(m) = std( ...
        per_repeat_table.NormalizedAUPRC(use_rows), ...
        0, ...
        'omitnan');

    mean_TopK(m) = mean( ...
        per_repeat_table.TopKRecovery(use_rows), ...
        'omitnan');

    sd_TopK(m) = std( ...
        per_repeat_table.TopKRecovery(use_rows), ...
        0, ...
        'omitnan');

    mean_BudgetFactor(m) = mean( ...
        per_repeat_table.Budget95Factor(use_rows), ...
        'omitnan');

    sd_BudgetFactor(m) = std( ...
        per_repeat_table.Budget95Factor(use_rows), ...
        0, ...
        'omitnan');

    mean_Precision95(m) = mean( ...
        per_repeat_table.PrecisionAt95Recall(use_rows), ...
        'omitnan');

    sd_Precision95(m) = std( ...
        per_repeat_table.PrecisionAt95Recall(use_rows), ...
        0, ...
        'omitnan');
end

summary_table = table( ...
    summary_method, ...
    mean_AUPRC, ...
    sd_AUPRC, ...
    mean_nAUPRC, ...
    sd_nAUPRC, ...
    mean_TopK, ...
    sd_TopK, ...
    mean_BudgetFactor, ...
    sd_BudgetFactor, ...
    mean_Precision95, ...
    sd_Precision95, ...
    'VariableNames', { ...
        'Method', ...
        'AUPRC_Mean', ...
        'AUPRC_SD', ...
        'NormalizedAUPRC_Mean', ...
        'NormalizedAUPRC_SD', ...
        'TopKRecovery_Mean', ...
        'TopKRecovery_SD', ...
        'Budget95Factor_Mean', ...
        'Budget95Factor_SD', ...
        'PrecisionAt95Recall_Mean', ...
        'PrecisionAt95Recall_SD'});

% Store candidate-list recovery at each budget.
candidate_rows = ...
    Nmethod * Nrep * numel(candidate_multipliers);

candidate_method     = strings(candidate_rows, 1);
candidate_repeat     = zeros(candidate_rows, 1);
candidate_multiplier = zeros(candidate_rows, 1);
candidate_recall_out = zeros(candidate_rows, 1);

row_id = 0;

for m = 1:Nmethod
    for rep = 1:Nrep
        for b = 1:numel(candidate_multipliers)

            row_id = row_id + 1;

            candidate_method(row_id) = ...
                methods(m).display_name;

            candidate_repeat(row_id) = rep;

            candidate_multiplier(row_id) = ...
                candidate_multipliers(b);

            candidate_recall_out(row_id) = ...
                candidate_recovery(b, m, rep);
        end
    end
end

candidate_table = table( ...
    candidate_method, ...
    candidate_repeat, ...
    candidate_multiplier, ...
    candidate_recall_out, ...
    'VariableNames', { ...
        'Method', ...
        'Repeat', ...
        'CandidateBudgetTimesK', ...
        'Recall'});

% Store mean precision values on the common recall grid.
PR_table = table( ...
    recall_grid, ...
    'VariableNames', {'Recall'});

for m = 1:Nmethod

    current_precision = squeeze( ...
        PR_precision_grid(:, m, :));

    mean_precision = mean( ...
        current_precision, ...
        2, ...
        'omitnan');

    variable_name = matlab.lang.makeValidName( ...
        methods(m).display_name);

    PR_table.(variable_name) = mean_precision;
end

%% Generate the two-panel comparison figure
method_labels = string({methods.display_name});

% Assign one color to each method.
plot_colors = lines(Nmethod);

% Draw HR-GGM last so overlapping curves remain visible.
plot_order = [2, 3, 4, 5, 6, 1];

line_styles = {'--', '-', '-.', ':', '-', '-'};
marker_styles = {'s', 'o', '^', 'd', 'v', 'p'};

% Set the final figure to 178 mm width and 1200 dpi.
FIGURE_DPI       = 1200;
FIGURE_WIDTH_IN  = 7.01;
FIGURE_HEIGHT_IN = 3.35;

AXIS_FONT_SIZE   = 6.3;
LABEL_FONT_SIZE  = 7.0;
TITLE_FONT_SIZE  = 7.3;
LEGEND_FONT_SIZE = 6.1;
REFERENCE_FONT_SIZE = 5.8;

line_widths = [1.90, 1.25, 1.25, 1.25, 1.25, 1.35];
marker_sizes = [4.4, 4.0, 4.0, 4.0, 4.0, 4.0];
ERROR_CAP_SIZE = 3;
REFERENCE_LINE_WIDTH = 0.85;

comparison_figure = figure( ...
    'Visible', 'off', ...
    'Color', 'w', ...
    'Units', 'inches', ...
    'Position', [0.5, 0.5, FIGURE_WIDTH_IN, FIGURE_HEIGHT_IN], ...
    'PaperUnits', 'inches', ...
    'PaperPosition', [0, 0, FIGURE_WIDTH_IN, FIGURE_HEIGHT_IN], ...
    'PaperSize', [FIGURE_WIDTH_IN, FIGURE_HEIGHT_IN], ...
    'Renderer', 'painters', ...
    'GraphicsSmoothing', 'on');

comparison_layout = tiledlayout( ...
    comparison_figure, ...
    1, ...
    2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

% Panel A: mean precision-recall curves.
ax1 = nexttile(comparison_layout, 1);
hold(ax1, 'on');

PR_handles = gobjects(Nmethod, 1);
mean_PR_curves = nan(numel(recall_grid), Nmethod);

for q = 1:numel(plot_order)

    m = plot_order(q);

    current_precision = squeeze( ...
        PR_precision_grid(:, m, :));

    mean_precision = mean( ...
        current_precision, ...
        2, ...
        'omitnan');

    mean_PR_curves(:, m) = mean_precision;

    PR_handles(m) = plot( ...
        ax1, ...
        recall_grid, ...
        mean_precision, ...
        'LineWidth', line_widths(m), ...
        'LineStyle', line_styles{m}, ...
        'Color', plot_colors(m, :), ...
        'DisplayName', methods(m).display_name);
end

xlim(ax1, [0, 1]);
ylim(ax1, [0, 1.015]);

xlabel(ax1, 'Recall', ...
    'FontSize', LABEL_FONT_SIZE);
ylabel(ax1, 'Precision', ...
    'FontSize', LABEL_FONT_SIZE);

title(ax1, ...
    '(A) Precision–recall', ...
    'FontSize', TITLE_FONT_SIZE, ...
    'FontWeight', 'bold');

grid(ax1, 'on');
box(ax1, 'off');

format_publication_axis(ax1, AXIS_FONT_SIZE);

hold(ax1, 'off');

% Panel B: candidate-list recovery with repeat-level variation.
ax2 = nexttile(comparison_layout, 2);
hold(ax2, 'on');

candidate_handles = gobjects(Nmethod, 1); %#ok<NASGU>

% Display-only x offsets can separate overlapping method markers.
x_offsets = ones(1, Nmethod);

for q = 1:numel(plot_order)

    m = plot_order(q);

    current_recovery = squeeze( ...
        candidate_recovery(:, m, :));

    current_mean = mean( ...
        current_recovery, ...
        2, ...
        'omitnan');

    current_sd = std( ...
        current_recovery, ...
        0, ...
        2, ...
        'omitnan');

    displayed_x = ...
        candidate_multipliers * x_offsets(m);

    candidate_handles(m) = errorbar( ...
        ax2, ...
        displayed_x, ...
        current_mean, ...
        current_sd, ...
        'LineStyle', line_styles{m}, ...
        'Marker', marker_styles{m}, ...
        'LineWidth', line_widths(m), ...
        'MarkerSize', marker_sizes(m), ...
        'CapSize', ERROR_CAP_SIZE, ...
        'Color', plot_colors(m, :), ...
        'MarkerFaceColor', 'none', ...
        'DisplayName', methods(m).display_name);
end

% Mark the candidate budget equal to the number of true edges.
hK = xline( ...
    ax2, ...
    1, ...
    '--', ...
    'K-edge budget', ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left', ...
    'LineWidth', REFERENCE_LINE_WIDTH, ...
    'HandleVisibility', 'off');
hK.FontName = 'Arial';
hK.FontSize = REFERENCE_FONT_SIZE;

% Mark the target recovery level.
hTarget = yline( ...
    ax2, ...
    target_recall, ...
    ':', ...
    sprintf('%.0f%% recovery', 100 * target_recall), ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'right', ...
    'LineWidth', REFERENCE_LINE_WIDTH, ...
    'HandleVisibility', 'off');
hTarget.FontName = 'Arial';
hTarget.FontSize = REFERENCE_FONT_SIZE;

set(ax2, 'XScale', 'log');

xlim(ax2, [ ...
    min(candidate_multipliers) * 0.85, ...
    max(candidate_multipliers) * 1.15]);

ylim(ax2, [0, 1.015]);

xticks(ax2, candidate_multipliers);
xticklabels(ax2, string(candidate_multipliers));

xlabel(ax2, 'Candidate-list size / K', ...
    'FontSize', LABEL_FONT_SIZE);
ylabel(ax2, 'Fraction of true edges recovered', ...
    'FontSize', LABEL_FONT_SIZE);

title(ax2, ...
    '(B) Candidate-list recovery', ...
    'FontSize', TITLE_FONT_SIZE, ...
    'FontWeight', 'bold');

grid(ax2, 'on');
box(ax2, 'off');

format_publication_axis(ax2, AXIS_FONT_SIZE);

hold(ax2, 'off');

% Use one shared legend for both panels.
comparison_legend = legend( ...
    ax1, ...
    PR_handles, ...
    method_labels, ...
    'NumColumns', 3, ...
    'Box', 'off', ...
    'FontName', 'Arial', ...
    'FontSize', LEGEND_FONT_SIZE);
comparison_legend.Layout.Tile = 'south';

drawnow;

%% Export figure files
figure_png = fullfile( ...
    output_dir, ...
    'Observed_variable_network_recovery.png');

figure_pdf = fullfile( ...
    output_dir, ...
    'Observed_variable_network_recovery.pdf');

figure_eps = fullfile( ...
    output_dir, ...
    'Observed_variable_network_recovery.eps');

figure_fig = fullfile( ...
    output_dir, ...
    'Observed_variable_network_recovery.fig');

exportgraphics( ...
    comparison_figure, ...
    figure_png, ...
    'Resolution', FIGURE_DPI, ...
    'BackgroundColor', 'white');

print( ...
    comparison_figure, ...
    figure_pdf, ...
    '-dpdf', ...
    '-painters');

print( ...
    comparison_figure, ...
    figure_eps, ...
    '-depsc', ...
    '-painters');

savefig(comparison_figure, figure_fig);
close(comparison_figure);

%% Save numerical results
excel_file = fullfile( ...
    output_dir, ...
    'Observed_variable_network_metrics.xlsx');

if exist(excel_file, 'file')
    delete(excel_file);
end

writetable( ...
    per_repeat_table, ...
    excel_file, ...
    'Sheet', 'PerRepeat');

writetable( ...
    summary_table, ...
    excel_file, ...
    'Sheet', 'Summary');

writetable( ...
    candidate_table, ...
    excel_file, ...
    'Sheet', 'CandidateRecovery');

writetable( ...
    PR_table, ...
    excel_file, ...
    'Sheet', 'MeanPRCurves');

%% Local functions
% Apply consistent axis formatting.
function format_publication_axis(ax, axis_font_size)

    ax.FontName = 'Arial';
    ax.FontSize = axis_font_size;
    ax.LineWidth = 0.7;
    ax.TickDir = 'out';
    ax.Box = 'off';
    ax.Layer = 'top';
    ax.GridAlpha = 0.18;
    ax.MinorGridAlpha = 0.10;
end

% Resolve one repeat-specific matrix file from supported extensions.
function matrix_file = resolve_matrix_file( ...
    folder_path, ...
    file_prefix, ...
    repeat_id)

    base_name = sprintf( ...
        '%s_rep%d', ...
        file_prefix, ...
        repeat_id);

    candidate_files = { ...
        fullfile(folder_path, [base_name, '.txt']), ...
        fullfile(folder_path, [base_name, '.csv']), ...
        fullfile(folder_path, [base_name, '.dat']), ...
        fullfile(folder_path, base_name)};

    matrix_file = '';

    for i = 1:numel(candidate_files)

        if exist(candidate_files{i}, 'file')
            matrix_file = candidate_files{i};
            return;
        end
    end

    possible_files = dir( ...
        fullfile(folder_path, [base_name, '.*']));

    possible_files = ...
        possible_files(~[possible_files.isdir]);

    if numel(possible_files) == 1

        matrix_file = fullfile( ...
            possible_files(1).folder, ...
            possible_files(1).name);

        return;

    elseif numel(possible_files) > 1

        possible_names = string( ...
            {possible_files.name});

        error(['Multiple files match %s in folder:\n%s\n\n', ...
               'Matching files:\n%s'], ...
               base_name, ...
               folder_path, ...
               strjoin(possible_names, newline));
    end

    error(['Cannot find matrix file %s in folder:\n%s'], ...
        base_name, ...
        folder_path);
end

% Read and validate a square numeric matrix.
function A = read_numeric_matrix(file_path)

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
        error('The matrix file is empty: %s', file_path);
    end

    valid_rows = ~all(isnan(A), 2);
    valid_cols = ~all(isnan(A), 1);

    A = A(valid_rows, valid_cols);

    if isempty(A)
        error('No numeric matrix was detected in: %s', file_path);
    end

    if size(A, 1) ~= size(A, 2)

        error(['The input matrix is not square after removing ', ...
               'nonnumeric rows and columns.\nFile: %s\n', ...
               'Detected size: %d x %d'], ...
               file_path, ...
               size(A, 1), ...
               size(A, 2));
    end

    A = double(A);
end

% Project a square matrix onto the symmetric matrix space.
function A = symmetrize_matrix(A)

    if size(A, 1) ~= size(A, 2)
        error('symmetrize_matrix requires a square matrix.');
    end

    A = (A + A') / 2;
end

% Calculate tie-aware precision-recall and candidate-budget metrics.
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

% Calculate the expected positives at a fixed budget under score ties.
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

% Calculate the expected budget required to reach a positive count.
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
