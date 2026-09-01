%% HR-GGM observed-network sparsity analysis
% Evaluates all undirected pairs of observed variables.
% The observed-variable count and simulation settings are read from the
% summary table. Metrics include baseline-corrected AUPRC, tie-aware
% top-K recovery, and relative off-diagonal coefficient error.
% Exported names use O, OH, and OG. Original input field names and
% source dataset identifiers are retained for compatibility and traceability.

clear;
clc;
close all force;

warning_state = warning;
warning('off', 'all');
warning_cleanup = onCleanup(@() warning(warning_state)); %#ok<NASGU>

script_file = mfilename('fullpath');
if isempty(script_file)
    root_dir = pwd;
else
    root_dir = fileparts(script_file);
end

%% Input and output paths

SUMMARY_FILE = fullfile(root_dir, 'All_conditions_summary.csv');
RESULT_ROOT = root_dir;
GT_ROOT = root_dir;

OUT_DIR = fullfile(RESULT_ROOT, 'O_Sparsity_Performance_Output');
if ~exist(OUT_DIR, 'dir')
    mkdir(OUT_DIR);
end

GT_TOL  = 1e-12;
EST_TOL = 1e-8;

%% Figure settings

STYLE.COMBINED_WIDTH_IN  = 7.20;
STYLE.COMBINED_HEIGHT_IN = 2.65;
STYLE.AXIS_FONT_SIZE   = 7.0;
STYLE.LABEL_FONT_SIZE  = 7.8;
STYLE.TITLE_FONT_SIZE  = 8.2;
STYLE.LEGEND_FONT_SIZE = 6.3;
STYLE.LINE_WIDTH       = 1.25;
STYLE.MARKER_SIZE      = 4.3;
STYLE.ERROR_CAP_SIZE   = 4;

STYLE.NOISE_COLORS = [
    0.0000 0.4470 0.7410
    0.8500 0.3250 0.0980
    0.9290 0.6940 0.1250
];

%% Dataset selection

T = readtable(SUMMARY_FILE, 'VariableNamingRule', 'preserve');

requiredVariables = { ...
    'Condition', 'ExperimentType', 'RepeatID', ...
    'Actual_LL_within_edges', 'Actual_LL_between_edges', ...
    'Actual_LL_total_edges', 'Noise_variance_fraction', ...
    'n_cells', 'n_local_genes', 'n_total_nodes'};

lls_require_variables(T, requiredVariables, SUMMARY_FILE);

T.Condition      = string(T.Condition);
T.ExperimentType = string(T.ExperimentType);

keep = T.ExperimentType == "LL_density" | ...
       T.ExperimentType == "Both_baseline";
Tll = T(keep, :);

% Calculate observed-network density and mean degree from the true edge count.
Tll.Plot_LL_total_edges = Tll.Actual_LL_total_edges;
Tll.LL_density_fraction = Tll.Actual_LL_total_edges ./ ...
    (Tll.n_local_genes .* (Tll.n_local_genes - 1) ./ 2);
Tll.Mean_LL_degree = 2 .* Tll.Actual_LL_total_edges ./ Tll.n_local_genes;

Tll = sortrows(Tll, ...
    {'Plot_LL_total_edges', 'Noise_variance_fraction', 'RepeatID'});

%% Dataset and design tables

listVars = { ...
    'Condition', 'ConditionLabel', 'ExperimentType', 'RepeatID', ...
    'n_cells', 'n_local_genes', 'n_total_nodes', ...
    'Actual_LL_within_edges', 'Actual_LL_between_edges', ...
    'Actual_LL_total_edges', 'Actual_LH_edges', 'Actual_LG_edges', ...
    'Noise_variance_fraction', 'LL_density_fraction', 'Mean_LL_degree'};
listVars = listVars(ismember(listVars, Tll.Properties.VariableNames));

DatasetList = Tll(:, listVars);
writetable(current_output_table(DatasetList), ...
    fullfile(OUT_DIR, 'O_sparsity_dataset_list_36.csv'));

DesignSummary = lls_make_design_summary(Tll);
writetable(current_output_table(DesignSummary), ...
    fullfile(OUT_DIR, 'O_sparsity_design_summary_12_conditions.csv'));

%% Dataset evaluation

metricRows = cell(0, 28);

for rr = 1:height(Tll)

    condition = Tll.Condition(rr);
    nLocal = Tll.n_local_genes(rr);
    nTotal = Tll.n_total_nodes(rr);

    resultFolder = lls_resolve_condition_folder(RESULT_ROOT, Tll, rr);
    gtFolder     = lls_resolve_condition_folder(GT_ROOT, Tll, rr);

    thetaEstFile = fullfile(resultFolder, 'build', 'Result', 'Theta_final.txt');
    thetaGTFile  = fullfile(gtFolder,     'build', 'Data',   'Theta_full.txt');

    if ~isfile(thetaEstFile)
        continue;
    end
    if ~isfile(thetaGTFile)
        continue;
    end

    ThetaEst = lls_read_theta(thetaEstFile, nTotal);
    ThetaGT  = lls_read_theta(thetaGTFile,  nTotal);

    if size(ThetaEst, 1) < nLocal || size(ThetaGT, 1) < nLocal
        continue;
    end

    LLest = ThetaEst(1:nLocal, 1:nLocal);
    LLgt  = ThetaGT(1:nLocal, 1:nLocal);

    % Evaluate all undirected observed-variable pairs together.
    mask = triu(true(nLocal), 1);
    estValues = LLest(mask);
    gtValues  = LLgt(mask);

    scores = abs(estValues);
    labels = abs(gtValues) > GT_TOL;

    N = numel(labels);
    K = sum(labels);

    auprc = lls_auprc_unique_thresholds(labels, scores);
    prevalence = K / N;
    normalizedAUPRC = (auprc - prevalence) / (1 - prevalence);

    [topKRecovery, expectedTP, boundaryTieCount, hasBoundaryTie] = ...
        lls_tie_aware_topK_recovery(labels, scores, K);

    % Relative Frobenius error of the complete off-diagonal observed block.
    LLestOff = LLest;
    LLgtOff  = LLgt;
    LLestOff(1:nLocal+1:end) = 0;
    LLgtOff(1:nLocal+1:end)  = 0;

    relativeLLError = norm(LLestOff - LLgtOff, 'fro') / ...
        max(norm(LLgtOff, 'fro'), eps);

    % Calculate coefficient and native-support diagnostics.
    trueEdgeEst = estValues(labels);
    trueEdgeGT  = gtValues(labels);

    if isempty(trueEdgeGT)
        trueEdgePearson = NaN;
        trueEdgeSignAccuracy = NaN;
    else
        if std(trueEdgeGT) > 0 && std(trueEdgeEst) > 0
            C = corrcoef(trueEdgeGT, trueEdgeEst);
            trueEdgePearson = C(1,2);
        else
            trueEdgePearson = NaN;
        end
        trueEdgeSignAccuracy = mean(sign(trueEdgeEst) == sign(trueEdgeGT));
    end

    nativeMask = scores > EST_TOL;
    nativeTP = sum(nativeMask & labels);
    nativeFP = sum(nativeMask & ~labels);
    nativeFN = sum(~nativeMask & labels);
    nativeEdgeCount = sum(nativeMask);

    nativePrecision = lls_safe_divide(nativeTP, nativeTP + nativeFP);
    nativeRecall    = lls_safe_divide(nativeTP, nativeTP + nativeFN);
    nativeF1        = lls_safe_divide(2 * nativePrecision * nativeRecall, ...
        nativePrecision + nativeRecall);

    metricRows(end+1, :) = { ...
        char(condition), ...
        char(lls_get_string(Tll, rr, 'ConditionLabel', "")), ...
        char(Tll.ExperimentType(rr)), ...
        Tll.RepeatID(rr), ...
        Tll.Noise_variance_fraction(rr), ...
        nLocal, ...
        nTotal, ...
        Tll.Actual_LL_within_edges(rr), ...
        Tll.Actual_LL_between_edges(rr), ...
        Tll.Actual_LL_total_edges(rr), ...
        Tll.LL_density_fraction(rr), ...
        Tll.Mean_LL_degree(rr), ...
        N, ...
        K, ...
        prevalence, ...
        auprc, ...
        normalizedAUPRC, ...
        topKRecovery, ...
        expectedTP, ...
        boundaryTieCount, ...
        hasBoundaryTie, ...
        relativeLLError, ...
        trueEdgePearson, ...
        trueEdgeSignAccuracy, ...
        nativeEdgeCount, ...
        nativePrecision, ...
        nativeRecall, ...
        nativeF1};
end

metricNames = { ...
    'Condition', 'ConditionLabel', 'ExperimentType', 'RepeatID', ...
    'Noise_variance_fraction', 'n_local_genes', 'n_total_nodes', ...
    'Actual_LL_within_edges', 'Actual_LL_between_edges', ...
    'Actual_LL_total_edges', 'LL_density_fraction', 'Mean_LL_degree', ...
    'N_candidate_LL_edges', 'N_true_LL_edges', 'Positive_fraction', ...
    'AUPRC', 'Normalized_AUPRC', 'TopK_recovery', ...
    'Expected_TP_at_K', 'Boundary_tie_count_at_K', ...
    'Has_boundary_tie_at_K', 'Relative_LL_Frobenius_error', ...
    'True_edge_Pearson_correlation', 'True_edge_sign_accuracy', ...
    'Native_LL_edge_count', 'Native_precision', ...
    'Native_recall', 'Native_F1'};

R = cell2table(metricRows, 'VariableNames', metricNames);

if isempty(R)
    error('No datasets were evaluated. Check RESULT_ROOT, GT_ROOT, and OutputFolder paths.');
end

writetable(current_output_table(R), fullfile(OUT_DIR, 'O_sparsity_metrics_each_dataset.csv'));

%% Summary across repeats

mainMetrics = { ...
    'AUPRC', ...
    'Normalized_AUPRC', ...
    'TopK_recovery', ...
    'Relative_LL_Frobenius_error', ...
    'True_edge_Pearson_correlation', ...
    'True_edge_sign_accuracy', ...
    'Native_LL_edge_count', ...
    'Native_precision', ...
    'Native_recall', ...
    'Native_F1'};

S = lls_group_mean_sd(R, ...
    {'Actual_LL_total_edges', 'Noise_variance_fraction'}, ...
    mainMetrics);

writetable(current_output_table(S), fullfile(OUT_DIR, 'O_sparsity_metrics_mean_sd.csv'));

%% Figure generation

lls_make_combined_figure(S, OUT_DIR, STYLE);

%% Helper functions

function lls_require_variables(T, names, filePath)
    missing = names(~ismember(names, T.Properties.VariableNames));
    if ~isempty(missing)
        error('Missing variables in %s: %s', filePath, strjoin(missing, ', '));
    end
end

function S = lls_make_design_summary(T)
    totalVals = sort(unique(T.Actual_LL_total_edges));
    noiseVals = sort(unique(T.Noise_variance_fraction));

    rows = cell(numel(totalVals) * numel(noiseVals), 13);
    r = 0;

    for ii = 1:numel(totalVals)
        for jj = 1:numel(noiseVals)
            idx = T.Actual_LL_total_edges == totalVals(ii) & ...
                  T.Noise_variance_fraction == noiseVals(jj);
            X = T(idx, :);
            if isempty(X)
                continue;
            end

            r = r + 1;
            rows(r,:) = { ...
                X.Actual_LL_within_edges(1), ...
                X.Actual_LL_between_edges(1), ...
                X.Actual_LL_total_edges(1), ...
                X.LL_density_fraction(1), ...
                X.Mean_LL_degree(1), ...
                X.Noise_variance_fraction(1), ...
                X.n_cells(1), ...
                X.n_local_genes(1), ...
                X.n_total_nodes(1), ...
                lls_first_numeric(X, 'Actual_LH_edges'), ...
                lls_first_numeric(X, 'Actual_LG_edges'), ...
                height(X), ...
                numel(unique(X.RepeatID))};
        end
    end

    rows = rows(1:r,:);
    names = { ...
        'Fixed_within_group_LL_edges', ...
        'Added_between_group_LL_edges', ...
        'Total_true_LL_edges', ...
        'LL_density_fraction', ...
        'Mean_LL_degree', ...
        'Noise_variance_fraction', ...
        'n_cells', ...
        'n_local_genes', ...
        'n_total_nodes', ...
        'True_LH_edges', ...
        'True_LG_edges', ...
        'Number_of_datasets', ...
        'Number_of_repeats'};

    S = cell2table(rows, 'VariableNames', names);
end

function x = lls_first_numeric(T, variableName)
    if ismember(variableName, T.Properties.VariableNames)
        x = T.(variableName)(1);
    else
        x = NaN;
    end
end

function folder = lls_resolve_condition_folder(rootDir, T, rowIndex)
    condition = char(string(T.Condition(rowIndex)));
    candidates = {fullfile(rootDir, condition)};

    if ismember('OutputFolder', T.Properties.VariableNames)
        raw = T.OutputFolder(rowIndex);
        if iscell(raw)
            raw = raw{1};
        end
        outputFolder = char(string(raw));
        if ~isempty(strtrim(outputFolder)) && ...
                ~strcmpi(strtrim(outputFolder), '<missing>')
            candidates{end+1} = outputFolder;
            candidates{end+1} = fullfile(rootDir, outputFolder);
        end
    end

    for ii = 1:numel(candidates)
        if isfolder(candidates{ii})
            folder = candidates{ii};
            return;
        end
    end

    folder = candidates{1};
end

function A = lls_read_theta(filePath, nExpected)
    X = readmatrix(filePath);
    if isempty(X)
        error('File is empty: %s', filePath);
    end

    X = X(~all(isnan(X), 2), :);

    if size(X,1) == size(X,2)
        A = X;
        A = lls_symmetrize(A);
        return;
    end

    if size(X,2) == 3
        ii = X(:,1);
        jj = X(:,2);
        vv = X(:,3);
        valid = ~isnan(ii) & ~isnan(jj) & ~isnan(vv);
        ii = ii(valid);
        jj = jj(valid);
        vv = vv(valid);

        if isempty(ii)
            A = zeros(nExpected, nExpected);
            return;
        end

        if min(ii) == 0 || min(jj) == 0
            ii = ii + 1;
            jj = jj + 1;
        end

        A = full(sparse(ii, jj, vv, nExpected, nExpected));
        A = lls_symmetrize(A);
        return;
    end

    error('Unrecognized matrix format: %s', filePath);
end

function A = lls_symmetrize(A)
    d = diag(A);
    A = 0.5 * (A + A');
    A(1:size(A,1)+1:end) = d;
end

function auprc = lls_auprc_unique_thresholds(labels, scores)
    labels = logical(labels(:));
    scores = scores(:);
    P = sum(labels);

    if P == 0
        auprc = NaN;
        return;
    end

    [scoresSorted, order] = sort(scores, 'descend');
    labelsSorted = labels(order);

    % Evaluate thresholds after complete tied-score groups.
    groupEnd = [find(diff(scoresSorted) ~= 0); numel(scoresSorted)];
    tpCum = cumsum(labelsSorted);

    tp = tpCum(groupEnd);
    selected = groupEnd;
    precision = tp ./ selected;
    recall = tp ./ P;

    % Add the no-edge point before trapezoidal integration of the PR curve.
    recall = [0; recall];
    precision = [1; precision];

    auprc = trapz(recall, precision);
end

function [recovery, expectedTP, tieCount, hasTie] = ...
        lls_tie_aware_topK_recovery(labels, scores, K)

    labels = logical(labels(:));
    scores = scores(:);
    P = sum(labels);

    if P == 0 || K <= 0
        recovery = NaN;
        expectedTP = NaN;
        tieCount = 0;
        hasTie = false;
        return;
    end

    K = min(K, numel(scores));
    sortedScores = sort(scores, 'descend');
    boundary = sortedScores(K);

    above = scores > boundary;
    tied  = scores == boundary;

    nAbove = sum(above);
    nNeedFromTie = K - nAbove;
    tieCount = sum(tied);
    hasTie = tieCount > nNeedFromTie;

    tpAbove = sum(labels & above);
    tpTied  = sum(labels & tied);

    expectedTP = tpAbove + nNeedFromTie * lls_safe_divide(tpTied, tieCount);
    recovery = expectedTP / P;
end

function value = lls_safe_divide(a, b)
    if isempty(b) || ~isfinite(b) || b == 0
        value = NaN;
    else
        value = a / b;
    end
end

function value = lls_get_string(T, rowIndex, variableName, defaultValue)
    if ismember(variableName, T.Properties.VariableNames)
        x = T.(variableName)(rowIndex);
        if iscell(x)
            x = x{1};
        end
        value = string(x);
    else
        value = string(defaultValue);
    end
end

function S = lls_group_mean_sd(T, groupVariables, metricVariables)
    keyTable = unique(T(:, groupVariables), 'rows', 'stable');
    nGroups = height(keyTable);

    S = keyTable;
    S.GroupCount = zeros(nGroups, 1);

    for mm = 1:numel(metricVariables)
        S.(['mean_' metricVariables{mm}]) = nan(nGroups, 1);
        S.(['std_' metricVariables{mm}])  = nan(nGroups, 1);
    end

    for gg = 1:nGroups
        idx = true(height(T),1);
        for kk = 1:numel(groupVariables)
            varName = groupVariables{kk};
            idx = idx & T.(varName) == keyTable.(varName)(gg);
        end

        S.GroupCount(gg) = sum(idx);

        for mm = 1:numel(metricVariables)
            varName = metricVariables{mm};
            x = T.(varName)(idx);
            x = x(isfinite(x));
            if isempty(x)
                mu = NaN;
                sd = NaN;
            else
                mu = mean(x);
                if numel(x) > 1
                    sd = std(x, 0);
                else
                    sd = 0;
                end
            end
            S.(['mean_' varName])(gg) = mu;
            S.(['std_' varName])(gg)  = sd;
        end
    end

    S = sortrows(S, groupVariables);
end

function lls_make_combined_figure(S, outDir, style)
    fig = figure('Visible','off', 'Color','w', 'Units','inches', ...
        'Position',[1 1 style.COMBINED_WIDTH_IN style.COMBINED_HEIGHT_IN], ...
        'PaperUnits','inches', ...
        'PaperPosition',[0 0 style.COMBINED_WIDTH_IN style.COMBINED_HEIGHT_IN], ...
        'PaperSize',[style.COMBINED_WIDTH_IN style.COMBINED_HEIGHT_IN]);

    tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

    nexttile;
    lls_plot_metric(S, 'Normalized_AUPRC', style, true);
    title('(A) Normalized AUPRC', 'FontSize',style.TITLE_FONT_SIZE, ...
        'FontWeight','bold');
    ylabel('Normalized AUPRC', 'FontSize',style.LABEL_FONT_SIZE);
    ylim(lls_probability_limits(S, 'Normalized_AUPRC'));

    nexttile;
    lls_plot_metric(S, 'TopK_recovery', style, false);
    title('(B) Top-\it{K} recovery', 'FontSize',style.TITLE_FONT_SIZE, ...
        'FontWeight','bold');
    ylabel('Observed-edge recovery', 'FontSize',style.LABEL_FONT_SIZE);
    ylim(lls_probability_limits(S, 'TopK_recovery'));

    nexttile;
    lls_plot_metric(S, 'Relative_LL_Frobenius_error', style, false);
    title('(C) Observed-block coefficient error', ...
        'FontSize',style.TITLE_FONT_SIZE, 'FontWeight','bold');
    ylabel('Relative Frobenius error', 'FontSize',style.LABEL_FONT_SIZE);
    ylim(lls_nonnegative_limits(S, 'Relative_LL_Frobenius_error'));

    lls_export(fig, outDir, 'Figure_O_sparsity_performance_1x3');
    close(fig);
end

function lls_plot_metric(S, metricName, style, showLegend)
    hold on;

    totalEdges = sort(unique(S.Actual_LL_total_edges));
    noiseVals  = sort(unique(S.Noise_variance_fraction));

    meanName = ['mean_' metricName];
    stdName  = ['std_' metricName];

    for ii = 1:numel(noiseVals)
        idx = S.Noise_variance_fraction == noiseVals(ii);
        X = S(idx,:);
        [~, order] = sort(X.Actual_LL_total_edges);
        X = X(order,:);

        x = zeros(height(X),1);
        for jj = 1:height(X)
            x(jj) = find(totalEdges == X.Actual_LL_total_edges(jj), 1);
        end

        colorIndex = min(ii, size(style.NOISE_COLORS,1));
        errorbar(x, X.(meanName), X.(stdName), '-o', ...
            'Color',style.NOISE_COLORS(colorIndex,:), ...
            'MarkerFaceColor',style.NOISE_COLORS(colorIndex,:), ...
            'LineWidth',style.LINE_WIDTH, ...
            'MarkerSize',style.MARKER_SIZE, ...
            'CapSize',style.ERROR_CAP_SIZE);
    end

    xlim([0.75 numel(totalEdges)+0.25]);
    xticks(1:numel(totalEdges));
    xticklabels(string(totalEdges));
    xlabel('True observed-network edges', ...
        'FontSize',style.LABEL_FONT_SIZE);

    grid on;
    box off;
    lls_format_axis(style);

    if showLegend
        labels = arrayfun(@(x) sprintf('Noise %.0f%%',100*x), ...
            noiseVals, 'UniformOutput',false);
        lgd = legend(labels, 'Location','best', ...
            'FontSize',style.LEGEND_FONT_SIZE);
        lgd.Box = 'off';
    end
end

function limits = lls_probability_limits(S, metricName)
    y = S.(['mean_' metricName]);
    e = S.(['std_' metricName]);
    lowValues = y - e;
    highValues = y + e;
    lowValues = lowValues(isfinite(lowValues));
    highValues = highValues(isfinite(highValues));

    if isempty(lowValues) || isempty(highValues)
        limits = [0 1];
        return;
    end

    low = min(lowValues);
    high = max(highValues);
    padding = max(0.01, 0.12 * max(high-low, 0.01));
    lower = max(0, low-padding);
    upper = min(1.005, max(1.0, high+padding));

    if upper-lower < 0.04
        lower = max(0, upper-0.04);
    end
    limits = [lower upper];
end

function limits = lls_nonnegative_limits(S, metricName)
    y = S.(['mean_' metricName]);
    e = S.(['std_' metricName]);
    highValues = y + e;
    highValues = highValues(isfinite(highValues));

    if isempty(highValues)
        high = NaN;
    else
        high = max(highValues);
    end

    if ~isfinite(high) || high <= 0
        limits = [0 1];
    else
        limits = [0 high*1.12];
    end
end

function lls_format_axis(style)
    set(gca, 'FontSize',style.AXIS_FONT_SIZE, ...
        'LineWidth',0.7, 'TickDir','out', ...
        'TickLength',[0.015 0.015]);
    ax = gca;
    ax.GridAlpha = 0.18;
end

function lls_export(fig, outDir, baseName)
    pngFile = fullfile(outDir, [baseName '.png']);
    pdfFile = fullfile(outDir, [baseName '.pdf']);
    figFile = fullfile(outDir, [baseName '.fig']);

    savefig(fig, figFile);
    try
        exportgraphics(fig, pngFile, 'Resolution',600);
        exportgraphics(fig, pdfFile, 'ContentType','vector');
    catch
        print(fig, pngFile, '-dpng', '-r600');
        print(fig, pdfFile, '-dpdf', '-painters');
    end
end

function out = current_output_table(T)
% Rename exported headers and descriptive labels without changing numeric data.
% Preserve source condition identifiers so they still match the input folders.

    out = T;
    names = out.Properties.VariableNames;
    names = regexprep(names, '(^|_)LL(?=_|$)', '$1O');
    names = regexprep(names, '(^|_)LH(?=_|$)', '$1OH');
    names = regexprep(names, '(^|_)LG(?=_|$)', '$1OG');
    names = regexprep(names, '^n_local_genes$', 'n_observed_variables');
    out.Properties.VariableNames = names;

    if ismember('ExperimentType', names)
        out.ExperimentType = replace(string(out.ExperimentType), ...
            "LL_density", "O_density");
    end
    if ismember('ConditionLabel', names)
        labels = string(out.ConditionLabel);
        labels = regexprep(labels, '(?<![A-Za-z])LL(?![A-Za-z])', 'O');
        labels = regexprep(labels, '(?<![A-Za-z])LH(?![A-Za-z])', 'OH');
        labels = regexprep(labels, '(?<![A-Za-z])LG(?![A-Za-z])', 'OG');
        out.ConditionLabel = labels;
    end
end
