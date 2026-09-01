%% Threshold-sensitive observed-network and hidden-regulator recovery
% Evaluates observed-network recovery with and without true hidden effects.
% Both settings are fitted with a hidden master regulator and a hidden
% global variable. Observed-regulatory recovery is evaluated only when
% true regulatory effects are present. Input folders are beside this script;
% outputs include per-dataset metrics, summaries, and a six-panel figure.
% Exported names use O and OH. Original dataset identifiers and paths
% are retained so existing input folders remain directly usable.

clear;
clc;
close all;

%% Figure defaults
% Apply a consistent publication font to all generated graphics.

set(groot, ...
    'defaultAxesFontName', 'Arial', ...
    'defaultTextFontName', 'Arial', ...
    'defaultLegendFontName', 'Arial');

%% Paths
% Use the script location as the root so the analysis remains portable.
script_file = mfilename('fullpath');

if isempty(script_file)
    ROOT_DIR = pwd;
else
    ROOT_DIR = fileparts(script_file);
end

OUT_DIR = fullfile( ...
    ROOT_DIR, ...
    'Sensitive_O_OH_Analysis');

%% Analysis settings
% Define the observed block O and hidden master regulatory variable H.
N_LOCAL  = 500;
N_HIDDEN = 1;

GT_TOL = 1e-12;

N_THRESHOLD_POINTS = 401;

LL_THRESHOLD_MAX_OVERRIDE = NaN;
LH_THRESHOLD_MAX_OVERRIDE = NaN;

THRESHOLD_PADDING = 1.05;

PLATEAU_FRACTION = 0.90;

FIGURE_DPI = 1200;

STYLE.FIG_WIDTH_IN  = 7.01;
STYLE.FIG_HEIGHT_IN = 4.80;

STYLE.AXIS_FONT_SIZE    = 6.3;
STYLE.LABEL_FONT_SIZE   = 7.0;
STYLE.TITLE_FONT_SIZE   = 7.3;
STYLE.LEGEND_FONT_SIZE  = 6.1;

STYLE.LINE_WIDTH        = 1.25;
STYLE.MARKER_SIZE       = 4.2;
STYLE.ERROR_CAP_SIZE    = 3;
STYLE.COLOR_FULL_HG = [0.0000, 0.4470, 0.7410];
STYLE.COLOR_LL_ONLY = [0.8500, 0.3250, 0.0980];

STYLE.COLOR_LH = [0.25, 0.25, 0.25];

%% Output preparation
% Create the analysis directory and replace the existing workbook.
if ~exist(ROOT_DIR, 'dir')
    error('ROOT_DIR does not exist: %s', ROOT_DIR);
end

if ~exist(OUT_DIR, 'dir')
    mkdir(OUT_DIR);
end

EXCEL_FILE = fullfile(OUT_DIR, 'Threshold_sensitive_results.xlsx');
if isfile(EXCEL_FILE)
    delete(EXCEL_FILE);
end

%% Dataset discovery
% Detect datasets with and without hidden effects and resolve their matrix files.
Manifest = scan_simulation_folders(ROOT_DIR);

if isempty(Manifest)
    error(['No matching datasets were found. Expected folder names ', ...
           'containing LLbetween or LLonly and noise_XXXpct_rep_XX.']);
end

Manifest = sort_manifest(Manifest);
writetable(Manifest, EXCEL_FILE, 'Sheet', 'Manifest');

validMask = Manifest.HasGroundTruth & Manifest.HasEstimate;
Manifest = Manifest(validMask, :);

if isempty(Manifest)
    error('No dataset contains both Theta_full.txt and Theta_final.txt.');
end

noiseLevelsAll = sort(unique(Manifest.NoisePct));
noiseColorsAll = make_noise_colors(numel(noiseLevelsAll));

%% Observed-network analysis
% Evaluate every undirected observed-variable pair on a common threshold grid.
LLRecords = load_LL_records(Manifest, N_LOCAL, GT_TOL);

if isempty(LLRecords)
    error('No observed-network dataset was loaded successfully.');
end

LL_THRESHOLD_MAX = choose_threshold_max( ...
    LLRecords, LL_THRESHOLD_MAX_OVERRIDE, THRESHOLD_PADDING, 'O');

LL_THRESHOLD_GRID = linspace(0, LL_THRESHOLD_MAX, N_THRESHOLD_POINTS)';

[LLMetrics, LLCurves] = analyze_records_by_threshold( ...
    LLRecords, LL_THRESHOLD_GRID, PLATEAU_FRACTION, false);

LLSummary = summarize_metric_table(LLMetrics, ...
    {'Setting', 'NoisePct'}, ...
    {'F1AtZero', 'MaximumF1', 'BestThreshold', ...
     'PrecisionAtBestThreshold', 'RecallAtBestThreshold', ...
     'AUPRC', 'NormalizedAUPRC', ...
     'PrecisionAtTrueEdgeBudget', ...
     'RecallAtTrueEdgeBudget', ...
     'F1AtTrueEdgeBudget', ...
     'BlockRelativeFrobeniusError'});

writetable(LLMetrics, fullfile(OUT_DIR, 'O_metrics_each_dataset.csv'));
writetable(LLCurves,  fullfile(OUT_DIR, 'O_threshold_curves_each_dataset.csv'));
writetable(LLSummary, fullfile(OUT_DIR, 'O_metrics_summary.csv'));

writetable(LLMetrics, EXCEL_FILE, 'Sheet', 'O_metrics');
writetable(LLCurves,  EXCEL_FILE, 'Sheet', 'O_curves');
writetable(LLSummary, EXCEL_FILE, 'Sheet', 'O_summary');

%% Observed-regulatory analysis
% Evaluate O-H entries only for datasets containing a true hidden master regulator.
ManifestFullHG = Manifest(Manifest.Setting == "Full_HG", :);

if isempty(ManifestFullHG)
    LHRecords = struct([]);
else
    LHRecords = load_LH_records(ManifestFullHG, N_LOCAL, N_HIDDEN, GT_TOL);

    if isempty(LHRecords)
        LHMetrics = table();
        LHCurves = table();
    else
        LH_THRESHOLD_MAX = choose_threshold_max( ...
            LHRecords, LH_THRESHOLD_MAX_OVERRIDE, THRESHOLD_PADDING, 'OH');

        LH_THRESHOLD_GRID = linspace(0, LH_THRESHOLD_MAX, ...
                                     N_THRESHOLD_POINTS)';

        [LHMetrics, LHCurves] = analyze_records_by_threshold( ...
            LHRecords, LH_THRESHOLD_GRID, PLATEAU_FRACTION, true);

        LHSummary = summarize_metric_table(LHMetrics, ...
            {'Setting', 'NoisePct'}, ...
            {'F1AtZero', 'MaximumF1', 'BestThreshold', ...
             'PrecisionAtBestThreshold', 'RecallAtBestThreshold', ...
             'AUPRC', 'NormalizedAUPRC', ...
             'PrecisionAtTrueEdgeBudget', ...
             'RecallAtTrueEdgeBudget', ...
             'F1AtTrueEdgeBudget', ...
             'SignInvariantRelativeFrobeniusError'});

        writetable(LHMetrics, ...
            fullfile(OUT_DIR, 'OH_metrics_each_dataset.csv'));
        writetable(LHCurves, ...
            fullfile(OUT_DIR, 'OH_threshold_curves_each_dataset.csv'));
        writetable(LHSummary, ...
            fullfile(OUT_DIR, 'OH_metrics_summary.csv'));

        writetable(LHMetrics, EXCEL_FILE, 'Sheet', 'OH_metrics');
        writetable(LHCurves,  EXCEL_FILE, 'Sheet', 'OH_curves');
        writetable(LHSummary, EXCEL_FILE, 'Sheet', 'OH_summary');

        noiseLevelsLH = sort(unique(LHMetrics.NoisePct));
        noiseColorsLH = colors_for_noise_subset( ...
            noiseLevelsAll, noiseColorsAll, noiseLevelsLH);

        make_combined_figure( ...
            LLCurves, LLMetrics, noiseLevelsAll, noiseColorsAll, ...
            LL_THRESHOLD_GRID, PLATEAU_FRACTION, ...
            LHCurves, LHMetrics, noiseLevelsLH, noiseColorsLH, ...
            LH_THRESHOLD_GRID, ...
            OUT_DIR, STYLE, FIGURE_DPI);
    end
end

%% Helper functions

function Manifest = scan_simulation_folders(rootDir)
% Identify condition folders and resolve ground-truth and estimated matrices.

    d = dir(rootDir);
    d = d([d.isdir]);
    d = d(~ismember({d.name}, {'.', '..'}));

    Condition       = strings(0, 1);
    Setting         = strings(0, 1);
    NoisePct        = zeros(0, 1);
    RepeatID        = zeros(0, 1);
    FolderPath      = strings(0, 1);
    GroundTruthPath = strings(0, 1);
    EstimatePath    = strings(0, 1);
    HasGroundTruth  = false(0, 1);
    HasEstimate     = false(0, 1);

    for i = 1:numel(d)
        name = string(d(i).name);
        low  = lower(name);

        if contains(low, 'llonly')
            setting = "Observed_only_HGmodel";
        elseif contains(low, 'llbetween') && contains(low, 'lh_')
            setting = "Full_HG";
        else
            continue;
        end

        token = regexpi(char(name), ...
            'noise_(\d+)pct_rep_(\d+)', 'tokens', 'once');

        if isempty(token)
            continue;
        end

        noisePct = str2double(token{1});
        repeatID = str2double(token{2});
        folder   = fullfile(rootDir, char(name));

        gtCandidates = { ...
            fullfile(folder, 'build', 'Data', 'Theta_full.txt'), ...
            fullfile(folder, 'Data', 'Theta_full.txt'), ...
            fullfile(folder, 'build', 'Theta_full.txt'), ...
            fullfile(folder, 'Theta_full.txt')};

        estCandidates = { ...
            fullfile(folder, 'build', 'Result', 'Theta_final.txt'), ...
            fullfile(folder, 'Result', 'Theta_final.txt'), ...
            fullfile(folder, 'build', 'Theta_final.txt'), ...
            fullfile(folder, 'Theta_final.txt')};

        [gtPath, hasGT]   = first_existing_file(gtCandidates);
        [estPath, hasEst] = first_existing_file(estCandidates);

        Condition(end+1, 1)       = name;
        Setting(end+1, 1)         = setting;
        NoisePct(end+1, 1)        = noisePct;
        RepeatID(end+1, 1)        = repeatID;
        FolderPath(end+1, 1)      = string(folder);
        GroundTruthPath(end+1, 1) = string(gtPath);
        EstimatePath(end+1, 1)    = string(estPath);
        HasGroundTruth(end+1, 1)  = hasGT;
        HasEstimate(end+1, 1)     = hasEst;
    end

    Manifest = table( ...
        Condition, Setting, NoisePct, RepeatID, FolderPath, ...
        GroundTruthPath, EstimatePath, HasGroundTruth, HasEstimate);
end

function Manifest = sort_manifest(Manifest)
% Order datasets by model setting, noise level, and repeat.

    order = 99 * ones(height(Manifest), 1);
    order(Manifest.Setting == "Full_HG")        = 1;
    order(Manifest.Setting == "Observed_only_HGmodel") = 2;

    Manifest.SettingOrder = order;
    Manifest = sortrows(Manifest, ...
        {'SettingOrder', 'NoisePct', 'RepeatID'});
    Manifest.SettingOrder = [];
end

function [filePath, existsFlag] = first_existing_file(candidates)
% Return the first existing path from a list of candidates.

    filePath = candidates{1};
    existsFlag = false;

    for i = 1:numel(candidates)
        if isfile(candidates{i})
            filePath = candidates{i};
            existsFlag = true;
            return;
        end
    end
end

function Records = load_LL_records(Manifest, nLocal, gtTol)
% Read and vectorize the symmetric observed blocks for all valid datasets.

    Records = struct( ...
        'Condition', {}, 'Setting', {}, 'NoisePct', {}, 'RepeatID', {}, ...
        'Labels', {}, 'Scores', {}, 'SignedGT', {}, 'SignedEst', {});

    for i = 1:height(Manifest)
        condition = Manifest.Condition(i);

        try
            ThetaGT  = read_theta_matrix( ...
                char(Manifest.GroundTruthPath(i)), []);
            ThetaEst = read_theta_matrix( ...
                char(Manifest.EstimatePath(i)), []);

            validate_LL_dimensions(ThetaGT, ThetaEst, nLocal, condition);

            LLgt  = symmetric_block(ThetaGT(1:nLocal, 1:nLocal));
            LLest = symmetric_block(ThetaEst(1:nLocal, 1:nLocal));

            [gtVec, estVec] = upper_triangle_vectors(LLgt, LLest);

            labels = abs(gtVec) > gtTol;
            scores = abs(estVec);

            if ~any(labels)
                continue;
            end

            r.Condition = condition;
            r.Setting   = Manifest.Setting(i);
            r.NoisePct  = Manifest.NoisePct(i);
            r.RepeatID  = Manifest.RepeatID(i);
            r.Labels    = labels(:);
            r.Scores    = scores(:);
            r.SignedGT  = gtVec(:);
            r.SignedEst = estVec(:);

            Records(end+1) = r; %#ok<AGROW>

        catch
            continue;
        end
    end
end

function Records = load_LH_records(Manifest, nLocal, nHidden, gtTol)
% Read and vectorize signed O-H blocks for datasets with hidden effects.

    Records = struct( ...
        'Condition', {}, 'Setting', {}, 'NoisePct', {}, 'RepeatID', {}, ...
        'Labels', {}, 'Scores', {}, 'SignedGT', {}, 'SignedEst', {});

    for i = 1:height(Manifest)
        condition = Manifest.Condition(i);

        try
            ThetaGT  = read_theta_matrix( ...
                char(Manifest.GroundTruthPath(i)), []);
            ThetaEst = read_theta_matrix( ...
                char(Manifest.EstimatePath(i)), []);

            validate_LH_dimensions( ...
                ThetaGT, ThetaEst, nLocal, nHidden, condition);

            gtVec  = signed_LH_vector(ThetaGT,  nLocal, nHidden);
            estVec = signed_LH_vector(ThetaEst, nLocal, nHidden);

            labels = abs(gtVec) > gtTol;
            scores = abs(estVec);

            if ~any(labels)
                continue;
            end

            r.Condition = condition;
            r.Setting   = Manifest.Setting(i);
            r.NoisePct  = Manifest.NoisePct(i);
            r.RepeatID  = Manifest.RepeatID(i);
            r.Labels    = labels(:);
            r.Scores    = scores(:);
            r.SignedGT  = gtVec(:);
            r.SignedEst = estVec(:);

            Records(end+1) = r; %#ok<AGROW>

        catch
            continue;
        end
    end
end

function thresholdMax = choose_threshold_max( ...
    Records, overrideValue, padding, blockName)
% Select a common upper bound for the absolute-threshold grid.

    if isfinite(overrideValue) && overrideValue > 0
        thresholdMax = overrideValue;
        return;
    end

    maxTrueScore = 0;
    maxAnyScore  = 0;

    for i = 1:numel(Records)
        labels = Records(i).Labels;
        scores = Records(i).Scores;

        if any(labels)
            maxTrueScore = max(maxTrueScore, max(scores(labels)));
        end

        maxAnyScore = max(maxAnyScore, max(scores));
    end

    base = maxTrueScore;

    if ~isfinite(base) || base <= 0
        base = maxAnyScore;
    end

    if ~isfinite(base) || base <= 0
        error('%s estimated scores are all zero or non-finite.', blockName);
    end

    thresholdMax = padding * base;
end

function [Metrics, Curves] = analyze_records_by_threshold( ...
    Records, thresholdGrid, plateauFraction, calculateLHDiagnostics)
% Calculate threshold curves and per-dataset recovery metrics.

    metricRows = cell(0, 30);
    curveRows  = cell(0, 1);

    for i = 1:numel(Records)
        R = Records(i);

        [Curve, M] = evaluate_threshold_curve( ...
            R.Labels, R.Scores, thresholdGrid, plateauFraction);

        blockTrueNorm = norm(R.SignedGT, 2);

        if blockTrueNorm > 0
            blockRelativeFrobeniusError = ...
                norm(R.SignedEst - R.SignedGT, 2) / blockTrueNorm;
        else
            blockRelativeFrobeniusError = NaN;
        end

        if calculateLHDiagnostics
            D = calculate_LH_diagnostics( ...
                R.Labels, R.SignedGT, R.SignedEst);
        else
            D = empty_LH_diagnostics();
        end

        metricRows(end+1, :) = { ... %#ok<AGROW>
            char(R.Condition), char(R.Setting), ...
            R.NoisePct, R.RepeatID, ...
            M.NCandidate, M.NTrue, M.PositiveFraction, ...
            M.F1AtZero, M.MaximumF1, M.BestThreshold, ...
            M.HighF1ThresholdLow, M.HighF1ThresholdHigh, ...
            M.HighF1ThresholdWidth, ...
            M.PrecisionAtBestThreshold, M.RecallAtBestThreshold, ...
            M.AUPRC, M.NormalizedAUPRC, ...
            M.PrecisionAtTrueEdgeBudget, ...
            M.RecallAtTrueEdgeBudget, ...
            M.F1AtTrueEdgeBudget, ...
            M.FPBeforeLastTrueEdge, ...
            M.LastTrueEdgeRankRatio, ...
            blockRelativeFrobeniusError, ...
            D.FalseSupportEnergyFraction, ...
            D.TrueSupportEnergyFraction, ...
            D.ScaleAlignedTrueSupportNRMSE, ...
            D.TrueSupportMagnitudeRatio, ...
            D.AbsolutePearsonTrueSupport, ...
            D.SignAgreementTrueSupport, ...
            D.SignInvariantRelativeFrobeniusError ...
            };

        nRows = height(Curve);

        newRows = table( ...
            repmat(string(R.Condition), nRows, 1), ...
            repmat(string(R.Setting), nRows, 1), ...
            repmat(R.NoisePct, nRows, 1), ...
            repmat(R.RepeatID, nRows, 1), ...
            Curve.Threshold, Curve.NPredicted, ...
            Curve.TP, Curve.FP, Curve.FN, ...
            Curve.Precision, Curve.Recall, Curve.F1, ...
            'VariableNames', { ...
                'Condition', 'Setting', 'NoisePct', 'RepeatID', ...
                'Threshold', 'NPredicted', 'TP', 'FP', 'FN', ...
                'Precision', 'Recall', 'F1'});

        curveRows{end+1, 1} = newRows; %#ok<AGROW>
    end

    metricNames = { ...
        'Condition', 'Setting', 'NoisePct', 'RepeatID', ...
        'NCandidate', 'NTrue', 'PositiveFraction', ...
        'F1AtZero', 'MaximumF1', 'BestThreshold', ...
        'HighF1ThresholdLow', 'HighF1ThresholdHigh', ...
        'HighF1ThresholdWidth', ...
        'PrecisionAtBestThreshold', 'RecallAtBestThreshold', ...
        'AUPRC', 'NormalizedAUPRC', ...
        'PrecisionAtTrueEdgeBudget', ...
        'RecallAtTrueEdgeBudget', ...
        'F1AtTrueEdgeBudget', ...
        'FPBeforeLastTrueEdge', ...
        'LastTrueEdgeRankRatio', ...
        'BlockRelativeFrobeniusError', ...
        'FalseSupportEnergyFraction', ...
        'TrueSupportEnergyFraction', ...
        'ScaleAlignedTrueSupportNRMSE', ...
        'TrueSupportMagnitudeRatio', ...
        'AbsolutePearsonTrueSupport', ...
        'SignAgreementTrueSupport', ...
        'SignInvariantRelativeFrobeniusError'};

    Metrics = cell2table(metricRows, 'VariableNames', metricNames);

    Metrics.Condition = string(Metrics.Condition);
    Metrics.Setting   = string(Metrics.Setting);

    for v = 3:numel(metricNames)
        name = metricNames{v};
        if iscell(Metrics.(name))
            Metrics.(name) = cell2mat(Metrics.(name));
        end
    end

    if isempty(curveRows)
        Curves = table();
    else
        Curves = vertcat(curveRows{:});
    end
end

function [Curve, M] = evaluate_threshold_curve( ...
    labels, scores, thresholdGrid, plateauFraction)
% Evaluate support recovery over one threshold grid.

    labels = logical(labels(:));
    scores = double(scores(:));
    scores(~isfinite(scores)) = 0;

    if numel(labels) ~= numel(scores)
        error('labels and scores must have equal length.');
    end

    nCandidate = numel(labels);
    nTrue      = sum(labels);

    if nTrue == 0
        error('The evaluated block contains no true edges.');
    end

    allScoresAsc  = sort(scores, 'ascend');
    trueScoresAsc = sort(scores(labels), 'ascend');

    nT = numel(thresholdGrid);

    nPredicted = zeros(nT, 1);
    TP         = zeros(nT, 1);
    FP         = zeros(nT, 1);
    FN         = zeros(nT, 1);
    precision  = nan(nT, 1);
    recall     = zeros(nT, 1);
    f1         = zeros(nT, 1);

    for j = 1:nT
        t = thresholdGrid(j);

        nPredicted(j) = count_strictly_greater_sorted(allScoresAsc, t);
        TP(j)         = count_strictly_greater_sorted(trueScoresAsc, t);
        FP(j)         = nPredicted(j) - TP(j);
        FN(j)         = nTrue - TP(j);

        if nPredicted(j) > 0
            precision(j) = TP(j) / nPredicted(j);
        end

        recall(j) = TP(j) / nTrue;

        if nPredicted(j) > 0 && ...
                (precision(j) + recall(j)) > 0
            f1(j) = 2 * precision(j) * recall(j) / ...
                    (precision(j) + recall(j));
        else
            f1(j) = 0;
        end
    end

    [maximumF1, maxIndices] = max(f1);
    maxIndices = find(abs(f1 - maximumF1) <= ...
        max(1e-12, 1e-10 * max(1, maximumF1)));

    bestIndex = maxIndices(round((numel(maxIndices) + 1) / 2));

    plateauMask = f1 >= plateauFraction * maximumF1;

    [leftIndex, rightIndex] = contiguous_interval_around_index( ...
        plateauMask, bestIndex);

    highLow   = thresholdGrid(leftIndex);
    highHigh  = thresholdGrid(rightIndex);
    highWidth = highHigh - highLow;

    prevalence = nTrue / nCandidate;
    auprc = threshold_based_auprc(labels, scores);

    if prevalence < 1
        normalizedAUPRC = (auprc - prevalence) / (1 - prevalence);
    else
        normalizedAUPRC = NaN;
    end

    [precisionAtK, recallAtK, f1AtK] = ...
        exact_true_edge_budget_metrics(labels, scores);

    weakestTrue = min(scores(labels));
    rankIncludingTies = sum(scores >= weakestTrue);
    fpBeforeLastTrue = rankIncludingTies - nTrue;
    lastTrueRankRatio = rankIncludingTies / nTrue;

    Curve = table( ...
        thresholdGrid(:), nPredicted, TP, FP, FN, ...
        precision, recall, f1, ...
        'VariableNames', { ...
            'Threshold', 'NPredicted', 'TP', 'FP', 'FN', ...
            'Precision', 'Recall', 'F1'});

    M.NCandidate = nCandidate;
    M.NTrue = nTrue;
    M.PositiveFraction = prevalence;
    M.F1AtZero = f1(1);
    M.MaximumF1 = maximumF1;
    M.BestThreshold = thresholdGrid(bestIndex);
    M.HighF1ThresholdLow = highLow;
    M.HighF1ThresholdHigh = highHigh;
    M.HighF1ThresholdWidth = highWidth;
    M.PrecisionAtBestThreshold = precision(bestIndex);
    M.RecallAtBestThreshold = recall(bestIndex);
    M.AUPRC = auprc;
    M.NormalizedAUPRC = normalizedAUPRC;
    M.PrecisionAtTrueEdgeBudget = precisionAtK;
    M.RecallAtTrueEdgeBudget = recallAtK;
    M.F1AtTrueEdgeBudget = f1AtK;
    M.FPBeforeLastTrueEdge = fpBeforeLastTrue;
    M.LastTrueEdgeRankRatio = lastTrueRankRatio;
end

function nGreater = count_strictly_greater_sorted(sortedAscending, threshold)
% Count sorted values that exceed a threshold.

    n = numel(sortedAscending);

    low  = 1;
    high = n;
    firstGreater = n + 1;

    while low <= high
        mid = floor((low + high) / 2);

        if sortedAscending(mid) > threshold
            firstGreater = mid;
            high = mid - 1;
        else
            low = mid + 1;
        end
    end

    nGreater = n - firstGreater + 1;
end

function [leftIndex, rightIndex] = contiguous_interval_around_index( ...
    logicalMask, centerIndex)
% Find the contiguous true interval containing a selected index.

    leftIndex = centerIndex;
    while leftIndex > 1 && logicalMask(leftIndex - 1)
        leftIndex = leftIndex - 1;
    end

    rightIndex = centerIndex;
    while rightIndex < numel(logicalMask) && logicalMask(rightIndex + 1)
        rightIndex = rightIndex + 1;
    end
end

function auprc = threshold_based_auprc(labels, scores)
% Sum precision weighted by recall increments after complete tied-score groups.

    labels = logical(labels(:));
    scores = scores(:);

    [scoresSorted, order] = sort(scores, 'descend');
    labelsSorted = labels(order);

    groupStart = [true; diff(scoresSorted) ~= 0];
    starts = find(groupStart);
    ends   = [starts(2:end) - 1; numel(scoresSorted)];

    groupSize = ends - starts + 1;
    groupTP   = zeros(numel(starts), 1);

    for g = 1:numel(starts)
        groupTP(g) = sum(labelsSorted(starts(g):ends(g)));
    end

    cumulativeTP = cumsum(groupTP);
    cumulativeN  = cumsum(groupSize);

    precision = cumulativeTP ./ cumulativeN;
    recall    = cumulativeTP ./ sum(labels);

    recallIncrement = diff([0; recall]);
    auprc = sum(recallIncrement .* precision);
end

function [precisionAtK, recallAtK, f1AtK] = ...
    exact_true_edge_budget_metrics(labels, scores)

    labels = logical(labels(:));
    scores = scores(:);

    K = sum(labels);

    [scoresSorted, order] = sort(scores, 'descend');
    labelsSorted = labels(order);

    groupStart = [true; diff(scoresSorted) ~= 0];
    starts = find(groupStart);
    ends   = [starts(2:end) - 1; numel(scoresSorted)];

    groupSize = ends - starts + 1;
    groupTP   = zeros(numel(starts), 1);

    for g = 1:numel(starts)
        groupTP(g) = sum(labelsSorted(starts(g):ends(g)));
    end

    cumulativeSize = cumsum(groupSize);
    cumulativeTP   = cumsum(groupTP);

    boundary = find(cumulativeSize >= K, 1, 'first');

    if boundary == 1
        previousSize = 0;
        previousTP   = 0;
    else
        previousSize = cumulativeSize(boundary - 1);
        previousTP   = cumulativeTP(boundary - 1);
    end

    slots = K - previousSize;
    positiveFraction = groupTP(boundary) / groupSize(boundary);

    expectedTP = previousTP + slots * positiveFraction;

    precisionAtK = expectedTP / K;
    recallAtK    = expectedTP / K;

    if precisionAtK + recallAtK > 0
        f1AtK = 2 * precisionAtK * recallAtK / ...
                (precisionAtK + recallAtK);
    else
        f1AtK = 0;
    end
end

function D = calculate_LH_diagnostics(labels, trueVector, estimatedVector)
% Calculate sign-invariant and support-specific O-H coefficient metrics.

    labels = logical(labels(:));
    trueVector = double(trueVector(:));
    estimatedVector = double(estimatedVector(:));

    trueSupport  = labels;
    falseSupport = ~labels;

    trueNorm = norm(trueVector, 2);

    if trueNorm > 0
        relativeErrorSameSign = ...
            norm(estimatedVector - trueVector, 2) / trueNorm;
        relativeErrorFlippedSign = ...
            norm(estimatedVector + trueVector, 2) / trueNorm;

        signInvariantRelativeFrobeniusError = ...
            min(relativeErrorSameSign, relativeErrorFlippedSign);
    else
        signInvariantRelativeFrobeniusError = NaN;
    end

    totalEnergy = sum(estimatedVector .^ 2);
    trueEnergy  = sum(estimatedVector(trueSupport) .^ 2);
    falseEnergy = sum(estimatedVector(falseSupport) .^ 2);

    if totalEnergy > 0
        falseFraction = falseEnergy / totalEnergy;
        trueFraction  = trueEnergy / totalEnergy;
    else
        falseFraction = NaN;
        trueFraction  = NaN;
    end

    gt = trueVector(trueSupport);
    est = estimatedVector(trueSupport);

    if norm(gt) > 0 && norm(est) > 0
        c = (est' * gt) / (est' * est);
        alignedEst = c * est;

        scaleAlignedNRMSE = norm(alignedEst - gt) / norm(gt);
        magnitudeRatio   = norm(est) / norm(gt);

        if numel(gt) >= 2 && std(gt) > 0 && std(est) > 0
            C = corrcoef(gt, est);
            absPearson = abs(C(1, 2));
        else
            absPearson = NaN;
        end

        globalSign = sign(gt' * est);
        if globalSign == 0
            globalSign = 1;
        end

        signAgreement = mean(sign(globalSign * est) == sign(gt));
    else
        scaleAlignedNRMSE = NaN;
        magnitudeRatio   = NaN;
        absPearson       = NaN;
        signAgreement    = NaN;
    end

    D.FalseSupportEnergyFraction = falseFraction;
    D.TrueSupportEnergyFraction = trueFraction;
    D.ScaleAlignedTrueSupportNRMSE = scaleAlignedNRMSE;
    D.TrueSupportMagnitudeRatio = magnitudeRatio;
    D.AbsolutePearsonTrueSupport = absPearson;
    D.SignAgreementTrueSupport = signAgreement;
    D.SignInvariantRelativeFrobeniusError = ...
        signInvariantRelativeFrobeniusError;
end

function D = empty_LH_diagnostics()
% Return missing-value placeholders for observed-network records.

    D.FalseSupportEnergyFraction = NaN;
    D.TrueSupportEnergyFraction = NaN;
    D.ScaleAlignedTrueSupportNRMSE = NaN;
    D.TrueSupportMagnitudeRatio = NaN;
    D.AbsolutePearsonTrueSupport = NaN;
    D.SignAgreementTrueSupport = NaN;
    D.SignInvariantRelativeFrobeniusError = NaN;
end

function Summary = summarize_metric_table(T, groupNames, metricNames)
% Aggregate selected metrics by setting and noise level.

    if isempty(T)
        Summary = table();
        return;
    end

    groupTable = unique(T(:, groupNames), 'rows', 'stable');
    rows = cell(height(groupTable), ...
        numel(groupNames) + 3 * numel(metricNames));

    outputNames = groupNames;

    for g = 1:height(groupTable)
        mask = true(height(T), 1);

        for j = 1:numel(groupNames)
            name = groupNames{j};

            if isstring(T.(name))
                mask = mask & T.(name) == groupTable.(name)(g);
            else
                mask = mask & T.(name) == groupTable.(name)(g);
            end

            rows{g, j} = groupTable.(name)(g);
        end

        col = numel(groupNames) + 1;

        for m = 1:numel(metricNames)
            metric = metricNames{m};
            values = T.(metric)(mask);
            values = values(isfinite(values));

            if isempty(values)
                meanValue = NaN;
                sdValue   = NaN;
                nValue    = 0;
            else
                meanValue = mean(values);
                sdValue   = std(values, 0);
                nValue    = numel(values);
            end

            rows{g, col}     = meanValue;
            rows{g, col + 1} = sdValue;
            rows{g, col + 2} = nValue;

            col = col + 3;
        end
    end

    for m = 1:numel(metricNames)
        metric = metricNames{m};

        outputNames{end+1} = ['Mean_' metric]; %#ok<AGROW>
        outputNames{end+1} = ['SD_' metric]; %#ok<AGROW>
        outputNames{end+1} = ['N_' metric]; %#ok<AGROW>
    end

    Summary = cell2table(rows, 'VariableNames', outputNames);

    for j = 1:numel(groupNames)
        name = groupNames{j};

        if isstring(T.(name))
            Summary.(name) = string(Summary.(name));
        elseif iscell(Summary.(name))
            Summary.(name) = cell2mat(Summary.(name));
        end
    end

    for j = (numel(groupNames) + 1):numel(outputNames)
        name = outputNames{j};
        if iscell(Summary.(name))
            Summary.(name) = cell2mat(Summary.(name));
        end
    end
end

function make_combined_figure( ...
    LLCurves, LLMetrics, noiseLevelsLL, noiseColorsLL, ...
    LLThresholdGrid, plateauFraction, ...
    LHCurves, LHMetrics, noiseLevelsLH, noiseColorsLH, ...
    LHThresholdGrid, outDir, style, dpi)
% Create and export the six-panel publication figure without displaying it.

    figWidth  = style.FIG_WIDTH_IN;
    figHeight = style.FIG_HEIGHT_IN;

    fig = figure( ...
        'Visible', 'off', ...
        'Color', 'w', ...
        'Units', 'inches', ...
        'Position', [0.8, 0.8, figWidth, figHeight], ...
        'PaperUnits', 'inches', ...
        'PaperPosition', [0, 0, figWidth, figHeight], ...
        'PaperSize', [figWidth, figHeight], ...
        'Renderer', 'painters', ...
        'GraphicsSmoothing', 'on');

    tl = tiledlayout(2, 3, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    LLxMax = effective_curve_axis_max( ...
        LLCurves, LLThresholdGrid, 0.01, 0.06);

    LHxMax = effective_curve_axis_max( ...
        LHCurves, LHThresholdGrid, 0.01, 0.05);

    % Panel A: observed-network recovery with regulatory and global effects.
    nexttile(1);

    plot_mean_threshold_curves( ...
        LLCurves, "Full_HG", ...
        noiseLevelsLL, noiseColorsLL, style);

    title('(A) With hidden effects', ...
        'FontSize', style.TITLE_FONT_SIZE, ...
        'FontWeight', 'bold');

    xlabel('Absolute truncation threshold', ...
        'FontSize', style.LABEL_FONT_SIZE);
    ylabel('Observed-network F1 score', ...
        'FontSize', style.LABEL_FONT_SIZE);

    xlim([0, LLxMax]);
    ylim([0, 1.02]);
    format_axis(style);

    legend(noise_legend_labels(noiseLevelsLL), ...
        'Location', 'northeast', ...
        'FontSize', style.LEGEND_FONT_SIZE, ...
        'Box', 'off');

    % Panel B: observed-network recovery without true hidden effects.
    nexttile(2);

    plot_mean_threshold_curves( ...
        LLCurves, "Observed_only_HGmodel", ...
        noiseLevelsLL, noiseColorsLL, style);

    title('(B) Without hidden effects', ...
        'FontSize', style.TITLE_FONT_SIZE, ...
        'FontWeight', 'bold');

    xlabel('Absolute truncation threshold', ...
        'FontSize', style.LABEL_FONT_SIZE);
    ylabel('Observed-network F1 score', ...
        'FontSize', style.LABEL_FONT_SIZE);

    xlim([0, LLxMax]);
    ylim([0, 1.02]);
    format_axis(style);

    legend(noise_legend_labels(noiseLevelsLL), ...
        'Location', 'northeast', ...
        'FontSize', style.LEGEND_FONT_SIZE, ...
        'Box', 'off');

    % Panel C: observed-network threshold intervals retaining at least 90% of maximum F1.
    nexttile(3);

    plot_LL_high_f1_intervals( ...
        LLMetrics, noiseLevelsLL, plateauFraction, style);

    title('(C) Observed-network thresholds', ...
        'FontSize', style.TITLE_FONT_SIZE, ...
        'FontWeight', 'bold');

    xlabel('Absolute truncation threshold', ...
        'FontSize', style.LABEL_FONT_SIZE);
    ylabel('Noise variance fraction', ...
        'FontSize', style.LABEL_FONT_SIZE);

    [xLow, xHigh] = interval_axis_limits( ...
        LLMetrics, LLThresholdGrid, 0.12);
    xlim([xLow, xHigh]);

    format_axis(style);

    % Panel D: recovery of observed-regulatory associations.
    nexttile(4);

    plot_mean_threshold_curves( ...
        LHCurves, "Full_HG", ...
        noiseLevelsLH, noiseColorsLH, style);

    title('(D) Observed-regulatory recovery', ...
        'FontSize', style.TITLE_FONT_SIZE, ...
        'FontWeight', 'bold');

    xlabel('Absolute truncation threshold', ...
        'FontSize', style.LABEL_FONT_SIZE);
    ylabel('Observed-regulatory F1 score', ...
        'FontSize', style.LABEL_FONT_SIZE);

    xlim([0, LHxMax]);
    ylim([0, 1.02]);
    format_axis(style);

    legend(noise_legend_labels(noiseLevelsLH), ...
        'Location', 'southwest', ...
        'FontSize', style.LEGEND_FONT_SIZE, ...
        'Box', 'off');

    % Panel E: off-diagonal relative coefficient error of Theta_O.
    nexttile(5);

    plot_metric_two_settings( ...
        LLMetrics, noiseLevelsLL, ...
        'BlockRelativeFrobeniusError', style);

    title('(E) Observed-block error', ...
        'FontSize', style.TITLE_FONT_SIZE, ...
        'FontWeight', 'bold');

    xlabel('Noise variance fraction', ...
        'FontSize', style.LABEL_FONT_SIZE);
    ylabel('Relative Frobenius error', ...
        'FontSize', style.LABEL_FONT_SIZE);

    format_nonnegative_error_axis( ...
        LLMetrics.BlockRelativeFrobeniusError);
    format_axis(style);

    legend({'With hidden effects', ...
            'Without hidden effects'}, ...
        'Location', 'best', ...
        'FontSize', style.LEGEND_FONT_SIZE, ...
        'Box', 'off');

    % Panel F: sign-invariant relative coefficient error of Theta_OH.
    nexttile(6);

    plot_metric_vs_noise( ...
        LHMetrics, noiseLevelsLH, ...
        'SignInvariantRelativeFrobeniusError', ...
        style.COLOR_LH, style);

    title('(F) Observed-regulatory error', ...
        'FontSize', style.TITLE_FONT_SIZE, ...
        'FontWeight', 'bold');

    xlabel('Noise variance fraction', ...
        'FontSize', style.LABEL_FONT_SIZE);
    ylabel('Relative Frobenius error', ...
        'FontSize', style.LABEL_FONT_SIZE);

    format_nonnegative_error_axis( ...
        LHMetrics.SignInvariantRelativeFrobeniusError);
    format_axis(style);

    drawnow;

    pngFile = fullfile(outDir, 'Figure_O_OH_combined_reordered_CleanLines.png');
    pdfFile = fullfile(outDir, 'Figure_O_OH_combined_reordered_CleanLines.pdf');
    figFile = fullfile(outDir, 'Figure_O_OH_combined_reordered_CleanLines.fig');
    savefig(fig, figFile);

    try
        exportgraphics(fig, pngFile, ...
            'Resolution', dpi, ...
            'BackgroundColor', 'white');
        exportgraphics(fig, pdfFile, ...
            'ContentType', 'vector', ...
            'BackgroundColor', 'white');
    catch
        print(fig, pngFile, '-dpng', sprintf('-r%d', dpi));
        print(fig, pdfFile, '-dpdf', '-painters');
    end

    close(fig);

end

function plot_LL_high_f1_intervals( ...
    Metrics, noiseLevels, plateauFraction, style)
% Plot repeat-level and mean high-F1 threshold intervals.

    hold on;

    yBase = 1:numel(noiseLevels);
    settingOffset = 0.16;
    repeatOffset  = [-0.045, 0, 0.045];

    settings   = ["Full_HG", "Observed_only_HGmodel"];
    labels     = {'With hidden effects', 'Without hidden effects'};
    colors     = [style.COLOR_FULL_HG; style.COLOR_LL_ONLY];
    directions = [-1, 1];

    for s = 1:numel(settings)
        setting = settings(s);

        for n = 1:numel(noiseLevels)
            noise = noiseLevels(n);

            rows = Metrics.Setting == setting & ...
                   Metrics.NoisePct == noise;

            T = Metrics(rows, :);
            T = sortrows(T, 'RepeatID');

            if isempty(T)
                continue;
            end

            yCenter = yBase(n) + directions(s) * settingOffset;

            for r = 1:height(T)
                if r <= numel(repeatOffset)
                    y = yCenter + repeatOffset(r);
                else
                    y = yCenter;
                end

                plot([T.HighF1ThresholdLow(r), T.HighF1ThresholdHigh(r)], ...
                     [y, y], ...
                     '-', ...
                     'Color', lighten_color(colors(s, :), 0.55), ...
                     'LineWidth', 1.0, ...
                     'HandleVisibility', 'off');

                plot(T.BestThreshold(r), y, 'o', ...
                    'MarkerSize', 4.0, ...
                    'MarkerFaceColor', 'w', ...
                    'MarkerEdgeColor', colors(s, :), ...
                    'LineWidth', 0.85, ...
                    'HandleVisibility', 'off');
            end

            meanLow  = mean(T.HighF1ThresholdLow,  'omitnan');
            meanHigh = mean(T.HighF1ThresholdHigh, 'omitnan');
            meanBest = mean(T.BestThreshold,       'omitnan');
            sdBest   = std( T.BestThreshold, 0,    'omitnan');

            plot([meanLow, meanHigh], [yCenter, yCenter], ...
                '-', 'Color', colors(s, :), ...
                'LineWidth', 3.0, ...
                'HandleVisibility', 'off');

            errorbar(meanBest, yCenter, 0, 0, sdBest, sdBest, ...
                'o', ...
                'Color', colors(s, :), ...
                'MarkerFaceColor', colors(s, :), ...
                'MarkerSize', style.MARKER_SIZE, ...
                'LineWidth', 1.2, ...
                'CapSize', style.ERROR_CAP_SIZE, ...
                'HandleVisibility', 'off');
        end
    end

    h1 = plot(nan, nan, '-', ...
        'Color', colors(1, :), 'LineWidth', 3.0);
    h2 = plot(nan, nan, '-', ...
        'Color', colors(2, :), 'LineWidth', 3.0);

    ax = gca;
    legendHandle = legend([h1, h2], labels, ...
        'Location', 'northeast', ...
        'FontSize', style.LEGEND_FONT_SIZE, ...
        'Box', 'off');

    yticks(yBase);
    yticklabels(compose('%g%%', noiseLevels));
    ylim([0.55, numel(noiseLevels) + 0.45]);

    % Place the legend in the gap between the two highest noise levels.
    % This avoids covering their threshold intervals while keeping the
    % legend inside Panel C.
    drawnow;
    ax.Units = 'normalized';
    legendHandle.Units = 'normalized';

    axisPosition = ax.Position;
    legendPosition = legendHandle.Position;
    yLimits = ylim(ax);
    gapCenter = numel(noiseLevels) - 0.5;
    gapFraction = (gapCenter - yLimits(1)) / diff(yLimits);

    legendPosition(1) = axisPosition(1) + axisPosition(3) - ...
        legendPosition(3) - 0.01 * axisPosition(3);
    legendPosition(2) = axisPosition(2) + ...
        gapFraction * axisPosition(4) - 0.5 * legendPosition(4);
    legendHandle.Position = legendPosition;

    grid on;
    box off;
end

function [xLow, xHigh] = interval_axis_limits( ...
    Metrics, thresholdGrid, paddingFraction)
% Set limits around the observed threshold intervals.

    values = [
        Metrics.HighF1ThresholdLow
        Metrics.HighF1ThresholdHigh
        Metrics.BestThreshold
    ];

    values = values(isfinite(values));

    gridLow  = thresholdGrid(1);
    gridHigh = thresholdGrid(end);

    if isempty(values)
        xLow  = gridLow;
        xHigh = gridHigh;
        return;
    end

    dataLow  = min(values);
    dataHigh = max(values);
    dataSpan = dataHigh - dataLow;

    if dataSpan <= 0
        dataSpan = max(0.05 * (gridHigh - gridLow), eps);
    end

    padding = paddingFraction * dataSpan;

    xLow  = max(gridLow,  dataLow  - padding);
    xHigh = min(gridHigh, dataHigh + padding);

    minimumWidth = max(0.08 * (gridHigh - gridLow), eps);

    if (xHigh - xLow) < minimumWidth
        center = 0.5 * (xLow + xHigh);
        xLow  = max(gridLow,  center - minimumWidth / 2);
        xHigh = min(gridHigh, center + minimumWidth / 2);

        if (xHigh - xLow) < minimumWidth
            if xLow <= gridLow
                xHigh = min(gridHigh, gridLow + minimumWidth);
            elseif xHigh >= gridHigh
                xLow = max(gridLow, gridHigh - minimumWidth);
            end
        end
    end
end

function plot_mean_threshold_curves( ...
    Curves, setting, noiseLevels, noiseColors, style)
% Plot mean F1 curves for each noise level.

    hold on;

    for n = 1:numel(noiseLevels)
        noise = noiseLevels(n);

        rows = Curves.Setting == setting & ...
               Curves.NoisePct == noise;

        if ~any(rows)
            continue;
        end

        T = Curves(rows, :);
        repeats   = sort(unique(T.RepeatID));
        thresholds = sort(unique(T.Threshold));

        Y = nan(numel(repeats), numel(thresholds));

        for r = 1:numel(repeats)
            Tr = T(T.RepeatID == repeats(r), :);
            Tr = sortrows(Tr, 'Threshold');

            [tf, loc] = ismember(thresholds, Tr.Threshold);
            Y(r, tf) = Tr.F1(loc(tf));
        end

        meanY = mean(Y, 1, 'omitnan');
        x = thresholds(:)';

        plot(x, meanY, ...
            'Color', noiseColors(n, :), ...
            'LineWidth', style.LINE_WIDTH, ...
            'LineStyle', '-');
    end

    grid on;
    box off;
end

function xMax = effective_curve_axis_max( ...
    Curves, thresholdGrid, minimumF1, paddingFraction)
% Set the displayed threshold range from non-negligible F1 values.

    relevant = Curves.Threshold( ...
        isfinite(Curves.F1) & Curves.F1 >= minimumF1);

    gridMax = thresholdGrid(end);

    if isempty(relevant)
        xMax = gridMax;
        return;
    end

    xMax = max(relevant) * (1 + paddingFraction);
    xMax = min(gridMax, xMax);

    minimumWidth = 0.20 * gridMax;
    xMax = max(xMax, minimumWidth);

    xMax = min(gridMax, nice_axis_upper(xMax));
end

function upper = nice_axis_upper(value)
% Round an axis maximum to a readable value.

    if ~isfinite(value) || value <= 0
        upper = 1;
        return;
    end

    magnitude = 10 ^ floor(log10(value));
    normalized = value / magnitude;

    if normalized <= 1
        nice = 1;
    elseif normalized <= 2
        nice = 2;
    elseif normalized <= 2.5
        nice = 2.5;
    elseif normalized <= 5
        nice = 5;
    else
        nice = 10;
    end

    upper = nice * magnitude;
end

function plot_metric_two_settings( ...
    Metrics, noiseLevels, metricName, style)
% Plot repeat values and mean with standard deviation for two settings.

    hold on;

    settings = ["Full_HG", "Observed_only_HGmodel"];
    colors = [style.COLOR_FULL_HG; style.COLOR_LL_ONLY];
    offsets = [-0.055, 0.055];

    for s = 1:numel(settings)
        meanValues = nan(numel(noiseLevels), 1);
        sdValues   = nan(numel(noiseLevels), 1);

        for n = 1:numel(noiseLevels)
            values = Metrics.(metricName)( ...
                Metrics.Setting == settings(s) & ...
                Metrics.NoisePct == noiseLevels(n));

            values = values(isfinite(values));

            if isempty(values)
                continue;
            end

            meanValues(n) = mean(values);
            sdValues(n)   = std(values, 0);

            jitter = linspace(-0.022, 0.022, numel(values));

            scatter( ...
                n + offsets(s) + jitter, values, ...
                28, ...
                'MarkerFaceColor', 'w', ...
                'MarkerEdgeColor', colors(s, :), ...
                'LineWidth', 0.9, ...
                'HandleVisibility', 'off');
        end

        errorbar( ...
            (1:numel(noiseLevels)) + offsets(s), ...
            meanValues, sdValues, ...
            '-o', ...
            'Color', colors(s, :), ...
            'MarkerFaceColor', colors(s, :), ...
            'MarkerSize', style.MARKER_SIZE, ...
            'LineWidth', style.LINE_WIDTH, ...
            'CapSize', style.ERROR_CAP_SIZE);
    end

    xlim([0.5, numel(noiseLevels) + 0.5]);
    xticks(1:numel(noiseLevels));
    xticklabels(compose('%g%%', noiseLevels));

    grid on;
    box off;
end

function format_nonnegative_error_axis(values)
% Set a nonnegative axis range for estimation errors.

    values = values(isfinite(values) & values >= 0);

    if isempty(values)
        ylim([0, 1]);
        return;
    end

    maximumValue = max(values);

    if maximumValue <= 0
        ylim([0, 0.1]);
        return;
    end

    upper = 1.15 * maximumValue;

    if upper <= 0.1
        step = 0.01;
    elseif upper <= 0.5
        step = 0.05;
    elseif upper <= 1
        step = 0.1;
    elseif upper <= 2
        step = 0.2;
    else
        magnitude = 10 ^ floor(log10(upper));
        step = magnitude / 2;
    end

    upper = ceil(upper / step) * step;
    ylim([0, upper]);
end

function plot_metric_vs_noise( ...
    Metrics, noiseLevels, metricName, color, style)
% Plot repeat values and mean with standard deviation across noise levels.

    hold on;

    meanValues = nan(numel(noiseLevels), 1);
    sdValues   = nan(numel(noiseLevels), 1);

    for n = 1:numel(noiseLevels)
        noise = noiseLevels(n);

        values = Metrics.(metricName)(Metrics.NoisePct == noise);
        values = values(isfinite(values));

        if isempty(values)
            continue;
        end

        meanValues(n) = mean(values);
        sdValues(n)   = std(values, 0);

        jitter = linspace(-0.055, 0.055, numel(values));

        scatter(n + jitter, values, ...
            28, ...
            'MarkerFaceColor', 'w', ...
            'MarkerEdgeColor', color, ...
            'LineWidth', 0.9, ...
            'HandleVisibility', 'off');
    end

    errorbar(1:numel(noiseLevels), ...
        meanValues, sdValues, ...
        '-o', ...
        'Color', color, ...
        'MarkerFaceColor', color, ...
        'MarkerSize', style.MARKER_SIZE, ...
        'LineWidth', style.LINE_WIDTH, ...
        'CapSize', style.ERROR_CAP_SIZE);

    xlim([0.5, numel(noiseLevels) + 0.5]);
    xticks(1:numel(noiseLevels));
    xticklabels(compose('%g%%', noiseLevels));

    grid on;
    box off;
end

function labels = noise_legend_labels(noiseLevels)
% Create noise-level legend labels.

    labels = cellstr(compose('%g%% noise', noiseLevels));
end

function format_axis(style)
% Apply publication formatting to the current axes.

    ax = gca;
    ax.FontName = 'Arial';
    ax.FontSize = style.AXIS_FONT_SIZE;
    ax.LineWidth = 0.7;
    ax.TickDir = 'out';
    ax.Box = 'off';
    ax.Layer = 'top';
end

function colors = make_noise_colors(n)
% Return distinct colors for the detected noise levels.

    base = [
        0.0000, 0.4470, 0.7410
        0.8500, 0.3250, 0.0980
        0.9290, 0.6940, 0.1250
        0.4940, 0.1840, 0.5560
        0.4660, 0.6740, 0.1880
        0.3010, 0.7450, 0.9330
        0.6350, 0.0780, 0.1840
    ];

    if n <= size(base, 1)
        colors = base(1:n, :);
    else
        colors = lines(n);
    end
end

function colorsSubset = colors_for_noise_subset( ...
    allNoise, allColors, subsetNoise)
% Map a subset of noise levels to the shared color scheme.

    colorsSubset = zeros(numel(subsetNoise), 3);

    for i = 1:numel(subsetNoise)
        idx = find(allNoise == subsetNoise(i), 1, 'first');

        if isempty(idx)
            colorsSubset(i, :) = [0.25, 0.25, 0.25];
        else
            colorsSubset(i, :) = allColors(idx, :);
        end
    end
end

function c = lighten_color(color, amount)
% Blend a color toward white.

    amount = max(0, min(1, amount));
    c = color + amount * (1 - color);
end

function B = symmetric_block(A)
% Average a square matrix with its transpose.

    B = 0.5 * (double(A) + double(A)');
end

function [gtVec, estVec] = upper_triangle_vectors(gtBlock, estBlock)
% Extract matched undirected upper-triangle entries.

    n = size(gtBlock, 1);
    mask = triu(true(n), 1);

    gtVec  = gtBlock(mask);
    estVec = estBlock(mask);
end

function vector = signed_LH_vector(Theta, nLocal, nHidden)
% Average O-H and transposed H-O entries while preserving sign.

    hiddenIdx = nLocal + (1:nHidden);

    blockLH = double(Theta(1:nLocal, hiddenIdx));
    blockHL = double(Theta(hiddenIdx, 1:nLocal))';

    vector = 0.5 * (blockLH + blockHL);
    vector = vector(:);
end

function validate_LL_dimensions(ThetaGT, ThetaEst, nLocal, condition)
% Verify that observed blocks are available in both matrices.

    if size(ThetaGT, 1) < nLocal || size(ThetaGT, 2) < nLocal
        error('GT matrix is smaller than %d x %d for %s.', ...
            nLocal, nLocal, condition);
    end

    if size(ThetaEst, 1) < nLocal || size(ThetaEst, 2) < nLocal
        error('Estimated matrix is smaller than %d x %d for %s.', ...
            nLocal, nLocal, condition);
    end
end

function validate_LH_dimensions( ...
    ThetaGT, ThetaEst, nLocal, nHidden, condition)
% Verify that observed-regulatory blocks are available in both matrices.

    required = nLocal + nHidden;

    if size(ThetaGT, 1) < required || size(ThetaGT, 2) < required
        error('GT matrix is smaller than %d x %d for %s.', ...
            required, required, condition);
    end

    if size(ThetaEst, 1) < required || size(ThetaEst, 2) < required
        error('Estimated matrix is smaller than %d x %d for %s.', ...
            required, required, condition);
    end
end

function A = read_theta_matrix(filePath, nExpected)
% Read dense, flattened, or sparse-triplet matrix input.

    if ~isfile(filePath)
        error('File does not exist: %s', filePath);
    end

    X = readmatrix(filePath);

    if isempty(X)
        error('File is empty: %s', filePath);
    end

    X = X(~all(isnan(X), 2), :);
    X = X(:, ~all(isnan(X), 1));

    if isempty(X)
        error('No numeric values were found in: %s', filePath);
    end

    if size(X, 1) == size(X, 2)
        A = double(X);
        return;
    end

    if isvector(X)
        v = X(:);
        v = v(isfinite(v));

        n = round(sqrt(numel(v)));

        if n * n == numel(v)
            A = reshape(v, n, n);
            return;
        end
    end

    if size(X, 2) == 3
        ii = X(:, 1);
        jj = X(:, 2);
        vv = X(:, 3);

        valid = isfinite(ii) & isfinite(jj) & isfinite(vv);

        ii = ii(valid);
        jj = jj(valid);
        vv = vv(valid);

        if isempty(ii)
            if nargin < 2 || isempty(nExpected)
                error(['Sparse triplet is empty and matrix size was ', ...
                       'not provided: %s'], filePath);
            end

            A = zeros(nExpected, nExpected);
            return;
        end

        if min(ii) == 0 || min(jj) == 0
            ii = ii + 1;
            jj = jj + 1;
        end

        if nargin < 2 || isempty(nExpected)
            nExpected = max([ii; jj]);
        end

        A = full(sparse(ii, jj, vv, nExpected, nExpected));
        A = symmetrize_sparse_input(A);
        return;
    end

    error('Unrecognized matrix format: %s', filePath);
end

function A = symmetrize_sparse_input(A)
% Symmetrize sparse input without duplicating one-sided entries.

    diagonalPart = diag(diag(A));
    upperPart = triu(A, 1);
    lowerPart = tril(A, -1);

    hasUpper = nnz(upperPart) > 0;
    hasLower = nnz(lowerPart) > 0;

    if hasUpper && ~hasLower
        A = diagonalPart + upperPart + upperPart';
    elseif hasLower && ~hasUpper
        A = diagonalPart + lowerPart + lowerPart';
    else
        A = 0.5 * (A + A');
    end
end
