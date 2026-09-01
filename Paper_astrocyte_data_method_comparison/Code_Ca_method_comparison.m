%% Ca2+-associated gene-ranking comparison: HR-GGM, LV-GGM, and prior correlation
% Compare non-prior gene rankings based on the HR-GGM observed--hidden
% regulatory precision block, LV-GGM latent components, and a
% prior-correlation baseline. The 42 prior genes are used to select the
% LV-GGM component and are excluded from evaluation.
clear;
clc;
close all;
rng(20260720, 'twister');

%% 1. Analysis settings
% Resolve input files from the current working directory.
precisionFile = resolveExistingFile({'S3_Precision_matrix.txt'});
geneListFile  = resolveExistingFile({'S1_whole_gene_list.xlsx'});
referenceFile = resolveExistingFile({'S5_Ca_associated_reference_label.xlsx'});
lvggmFile = resolveExistingFile({'LVGGM_L_lowrank.txt'});
expressionFile = resolveExistingFile({'Data_whole.txt'});

% Workbook sheet names.
geneListSheet = 'Whole_genes';
seedGeneSheet = 'Ca_seed_genes';
referenceSheet = 'Annotation';

% Matrix dimensions and variable indices.
nGene = 3000;
nPriorSeeds = 42;
idxHiddenCa = 3001;
idxGlobal = 3002; %#ok<NASGU>
expectedPrecisionSize = nGene + 2;

% Ranking budgets.
KList = [10, 20, 30, 50, 75, 100, 150, 200, 300, ...
         500, 721, 750, 1000, nGene - nPriorSeeds].';
fixedTopK = 721;

% Spectral score-extraction settings.
relativeEigenTolerance = 1e-8;
absoluteEigenTolerance = 1e-12;
seedDirectionAlpha = 0.05;
nSeedPermutations = 500;
maxSeedDirections = 20;
seedPcaEnergyTarget = 0.90;

% Resampling settings for external-label inference.
nLabelPermutations = 5000;
nBootstrap = 3000;
bootstrapAlpha = 0.05;

% Optional label-informed sensitivity diagnostic.
runOracleDiagnostic = true;
oracleMaxComponents = 100;

% Output directory.
outputFolder = 'Method_Comparison_Output';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

fprintf('============================================================\n');
fprintf('HR-GGM, LV-GGM, and prior-correlation comparison\n');
fprintf('============================================================\n');
fprintf('Output folder: %s\n\n', outputFolder);
fprintf('LV-GGM latent-dimension whitening: OFF\n\n');

%% 2. Gene names, prior genes, and reference labels
fprintf('Reading gene-name workbook...\n');
geneOpts = detectImportOptions(geneListFile, 'Sheet', geneListSheet, ...
    'VariableNamingRule', 'preserve');
GeneList = readtable(geneListFile, geneOpts);

symbolVar = findVariableName(GeneList, ...
    {'SYMBOL', 'Symbol', 'Gene', 'Gene Symbol'});
geneIDVar = findVariableName(GeneList, ...
    {'ENSMUSG ID', 'ENSMUSG_ID', 'Ensembl ID', 'EnsemblID'});

matrixGeneNames = strtrim(string(GeneList.(symbolVar)));
matrixGeneIDs = strtrim(string(GeneList.(geneIDVar)));

if height(GeneList) ~= nGene
    error('Whole_genes must contain exactly %d genes; found %d.', ...
        nGene, height(GeneList));
end
if any(strlength(matrixGeneNames) == 0) || any(strlength(matrixGeneIDs) == 0)
    error('Whole_genes contains empty gene symbols or Ensembl IDs.');
end
if numel(unique(matrixGeneNames)) ~= nGene
    error('Whole_genes contains duplicate gene symbols.');
end
if numel(unique(matrixGeneIDs)) ~= nGene
    error('Whole_genes contains duplicate Ensembl IDs.');
end

seedOpts = detectImportOptions(geneListFile, 'Sheet', seedGeneSheet, ...
    'VariableNamingRule', 'preserve');
SeedList = readtable(geneListFile, seedOpts);
seedSymbolVar = findVariableName(SeedList, ...
    {'SYMBOL', 'Symbol', 'Gene', 'Gene Symbol'});
seedGeneNames = strtrim(string(SeedList.(seedSymbolVar)));

if height(SeedList) ~= nPriorSeeds
    error('Ca_seed_genes must contain exactly %d genes; found %d.', ...
        nPriorSeeds, height(SeedList));
end
if ~isequal(sort(matrixGeneNames(1:nPriorSeeds)), sort(seedGeneNames))
    error(['The first %d matrix genes are not the same set as the ' ...
           'Ca_seed_genes sheet.'], nPriorSeeds);
end

fprintf('Reading independent Mark-T/F annotation...\n');
refOpts = detectImportOptions(referenceFile, 'Sheet', referenceSheet, ...
    'VariableNamingRule', 'preserve');
Ref = readtable(referenceFile, refOpts);

geneVar = findVariableName(Ref, {'Gene'});
markVar = findVariableName(Ref, ...
    {'Mark T/F', 'MarkT/F', 'Mark_T_F', 'MarkTF', 'Mark T F'});
weightVar = findVariableName(Ref, ...
    {'AbsWeight', 'Abs Weight', 'AbsoluteWeight', 'Absolute Weight'});

referenceGeneNames = strtrim(string(Ref.(geneVar)));
markLabel = upper(strtrim(string(Ref.(markVar))));
referenceAbsWeight = numericColumn(Ref.(weightVar));
validRows = strlength(referenceGeneNames) > 0 & ...
            (markLabel == "T" | markLabel == "F") & ...
            isfinite(referenceAbsWeight);
if any(~validRows)
    fprintf('Removing %d invalid reference rows.\n', sum(~validRows));
    Ref = Ref(validRows, :);
    referenceGeneNames = referenceGeneNames(validRows);
    markLabel = markLabel(validRows);
    referenceAbsWeight = referenceAbsWeight(validRows);
end

nEval = nGene - nPriorSeeds;
if height(Ref) ~= nEval
    error('Annotation must contain exactly %d non-prior genes; found %d.', ...
        nEval, height(Ref));
end

expectedEvalNames = matrixGeneNames(nPriorSeeds + 1:nGene);
expectedEvalIDs = matrixGeneIDs(nPriorSeeds + 1:nGene);
nameMatch = expectedEvalNames == referenceGeneNames;
if any(~nameMatch)
    j = find(~nameMatch, 1, 'first');
    error(['Gene-order mismatch at annotation row %d / matrix row %d: ' ...
           'expected "%s", found "%s".'], ...
          j, j + nPriorSeeds, expectedEvalNames(j), referenceGeneNames(j));
end

isMarkT = markLabel == "T";
isMarkF = markLabel == "F";
nMarkT = sum(isMarkT);
nMarkF = sum(isMarkF);
prevalence = nMarkT / nEval;
originalEvalIndex = (1:nEval).';
matrixEvalIndex = (nPriorSeeds + 1:nGene).';

fprintf('  Gene order verified for all %d evaluation genes.\n', nEval);
fprintf('  Mark-T = %d; Mark-F = %d; prevalence = %.6f.\n\n', ...
    nMarkT, nMarkF, prevalence);

%% 3. HR-GGM Ca2+-associated scores
fprintf('Reading HR-GGM precision matrix...\n');
Theta = readmatrix(precisionFile);
if ~isequal(size(Theta), [expectedPrecisionSize, expectedPrecisionSize])
    error('HR-GGM precision matrix must be %d x %d; found %d x %d.', ...
        expectedPrecisionSize, expectedPrecisionSize, ...
        size(Theta, 1), size(Theta, 2));
end
if any(~isfinite(Theta), 'all')
    error('HR-GGM precision matrix contains NaN or Inf.');
end

hrSymmetryError = max(abs(Theta - Theta.'), [], 'all');
allHrScores = abs(Theta(1:nGene, idxHiddenCa));
hrScore = allHrScores(nPriorSeeds + 1:nGene);
clear Theta allHrScores;

maxHrWeightDifference = max(abs(hrScore - referenceAbsWeight));
if maxHrWeightDifference > 1e-10
    error(['The HR-GGM matrix-derived Ca weights do not match the ' ...
           'AbsWeight column. Maximum difference = %.3e.'], ...
          maxHrWeightDifference);
end

fprintf('  HR-GGM symmetry error: %.3e\n', hrSymmetryError);
fprintf('  Matrix-versus-Excel AbsWeight difference: %.3e\n', ...
    maxHrWeightDifference);
fprintf('  Extracted |Theta(gene, hidden Ca)| for %d non-prior genes.\n\n', nEval);

%% 4. Latent-matrix score extraction
% Parameters passed to the LV-GGM latent-matrix analysis.
latentOptions = struct();
latentOptions.nGene = nGene;
latentOptions.seedIndex = (1:nPriorSeeds).';
latentOptions.evalIndex = (nPriorSeeds + 1:nGene).';
latentOptions.relativeEigenTolerance = relativeEigenTolerance;
latentOptions.absoluteEigenTolerance = absoluteEigenTolerance;
latentOptions.seedDirectionAlpha = seedDirectionAlpha;
latentOptions.nSeedPermutations = nSeedPermutations;
latentOptions.maxSeedDirections = maxSeedDirections;
latentOptions.seedPcaEnergyTarget = seedPcaEnergyTarget;
latentOptions.runOracleDiagnostic = runOracleDiagnostic;
latentOptions.oracleMaxComponents = oracleMaxComponents;
latentOptions.evalLabels = isMarkT;
latentOptions.originalEvalIndex = originalEvalIndex;

fprintf('Processing LV-GGM latent matrix...\n');
[lvScores, lvDiag, lvSpectrum, lvDirectionNull] = ...
    processLatentMatrix(lvggmFile, 'LV-GGM', latentOptions);

fprintf('\nComputing prior-correlation baseline...\n');
priorCorrScore = computePriorCorrelationScore(expressionFile, nGene, nPriorSeeds);

% Methods included in all numerical and graphical comparisons.
methodNames = ["HR-GGM"; "LV-GGM"; "Prior correlation"];
methodScores = [hrScore, ...
                lvScores.PriorSelectedComponent, ...
                priorCorrScore];
nMethod = numel(methodNames);

figureMethodNames = methodNames;
figureMethodScores = methodScores;
nFigureMethod = nMethod;

%% 5. Complete-ranking evaluation
MetricRows = cell(nMethod, 1);
TopKLong = table();
RankedLong = table();
ContingencyRows = cell(nMethod, 1);
PR = cell(nMethod, 1);

for m = 1:nMethod
    evalResult = evaluateRanking(methodScores(:, m), isMarkT, KList, ...
        fixedTopK, originalEvalIndex);

    MetricRows{m} = table(methodNames(m), evalResult.AP, ...
        evalResult.NormalizedAP, evalResult.AUROC, ...
        evalResult.PositiveRateAtFixedK, ...
        evalResult.EnrichmentAtFixedK, ...
        evalResult.RecallAtFixedK, evalResult.OddsRatio, ...
        evalResult.OR_CI95(1), evalResult.OR_CI95(2), ...
        evalResult.FisherRightP, evalResult.CutoffScore, ...
        evalResult.NextScore, evalResult.TiedAtCutoff, ...
        'VariableNames', {'Method', 'AveragePrecision', ...
        'NormalizedAP', 'AUROC', 'PositiveRateTop721', ...
        'EnrichmentTop721', 'RecallTop721', 'OddsRatioTop721', ...
        'OR_CI95_Lower', 'OR_CI95_Upper', 'FisherRightP', ...
        'CutoffScoreTop721', 'NextScore', 'TiedAtCutoff'});

    tmpTopK = table(repmat(methodNames(m), numel(KList), 1), KList, ...
        evalResult.TopKPositiveCount, evalResult.TopKPositiveRate, ...
        evalResult.TopKRecall, evalResult.TopKEnrichment, ...
        'VariableNames', {'Method', 'TopK', 'MarkTCount', ...
        'PositiveRate', 'Recall', 'EnrichmentFold'});
    TopKLong = [TopKLong; tmpTopK]; %#ok<AGROW>

    rankN = (1:nEval).';
    tmpRank = table(repmat(methodNames(m), nEval, 1), rankN, ...
        matrixEvalIndex(evalResult.RankOrder), ...
        expectedEvalIDs(evalResult.RankOrder), ...
        referenceGeneNames(evalResult.RankOrder), ...
        markLabel(evalResult.RankOrder), ...
        methodScores(evalResult.RankOrder, m), ...
        cumsum(isMarkT(evalResult.RankOrder)), ...
        evalResult.RunningPositiveRate, evalResult.RunningEnrichment, ...
        'VariableNames', {'Method', 'Rank', 'MatrixIndex', 'ENSMUSG_ID', ...
        'Gene', 'ReferenceLabel', 'Score', 'CumulativeMarkT', ...
        'RunningPositiveRate', 'RunningEnrichment'});
    RankedLong = [RankedLong; tmpRank]; %#ok<AGROW>

    ContingencyRows{m} = table(methodNames(m), ...
        evalResult.Contingency(1,1), evalResult.Contingency(1,2), ...
        evalResult.Contingency(2,1), evalResult.Contingency(2,2), ...
        evalResult.OddsRatio, evalResult.OR_CI95(1), ...
        evalResult.OR_CI95(2), evalResult.FisherRightP, ...
        'VariableNames', {'Method', 'MarkT_ModelT', 'MarkT_ModelF', ...
        'MarkF_ModelT', 'MarkF_ModelF', 'OddsRatio', ...
        'OR_CI95_Lower', 'OR_CI95_Upper', 'FisherRightP'});

    PR{m} = evalResult.PR;
end

MethodMetrics = vertcat(MetricRows{:});
ContingencyStats = vertcat(ContingencyRows{:});

%% 6. Plotted-method evaluation
FigureMetricRows = cell(nFigureMethod, 1);
FigureTopKLong = table();
FigurePR = cell(nFigureMethod, 1);

for m = 1:nFigureMethod
    figEvalResult = evaluateRanking(figureMethodScores(:, m), isMarkT, KList, ...
        fixedTopK, originalEvalIndex);

    FigureMetricRows{m} = table(figureMethodNames(m), figEvalResult.AP, ...
        figEvalResult.NormalizedAP, figEvalResult.AUROC, ...
        figEvalResult.PositiveRateAtFixedK, ...
        figEvalResult.EnrichmentAtFixedK, ...
        figEvalResult.RecallAtFixedK, figEvalResult.OddsRatio, ...
        figEvalResult.OR_CI95(1), figEvalResult.OR_CI95(2), ...
        figEvalResult.FisherRightP, ...
        'VariableNames', {'Method', 'AveragePrecision', ...
        'NormalizedAP', 'AUROC', 'PositiveRateTop721', ...
        'EnrichmentTop721', 'RecallTop721', 'OddsRatioTop721', ...
        'OR_CI95_Lower', 'OR_CI95_Upper', 'FisherRightP'});

    tmpTopK = table(repmat(figureMethodNames(m), numel(KList), 1), KList, ...
        figEvalResult.TopKPositiveCount, figEvalResult.TopKPositiveRate, ...
        figEvalResult.TopKRecall, figEvalResult.TopKEnrichment, ...
        'VariableNames', {'Method', 'TopK', 'MarkTCount', ...
        'PositiveRate', 'Recall', 'EnrichmentFold'});
    FigureTopKLong = [FigureTopKLong; tmpTopK]; %#ok<AGROW>

    FigurePR{m} = figEvalResult.PR;
end

FigureMethodMetrics = vertcat(FigureMetricRows{:});

%% 7. Label-permutation tests
fprintf('\nRunning %d shared label permutations...\n', nLabelPermutations);
permAP = zeros(nLabelPermutations, nMethod);
permEnrich721 = zeros(nLabelPermutations, nMethod);
observedAP = MethodMetrics.AveragePrecision.';
observedEnrich721 = MethodMetrics.EnrichmentTop721.';

% Fixed method-specific ranking orders for shared label permutations.
rankOrders = zeros(nEval, nMethod);
for m = 1:nMethod
    rankOrders(:,m) = deterministicRankOrder(methodScores(:,m), originalEvalIndex);
end

for b = 1:nLabelPermutations
    yPerm = isMarkT(randperm(nEval));
    yPerm = yPerm(:);
    for m = 1:nMethod
        yRank = yPerm(rankOrders(:,m));
        permAP(b,m) = averagePrecisionFromRankedLabels(yRank);
        permEnrich721(b,m) = mean(yRank(1:fixedTopK)) / prevalence;
    end
end

apPermP = (1 + sum(permAP >= observedAP, 1)) ./ (nLabelPermutations + 1);
enrichPermP = (1 + sum(permEnrich721 >= observedEnrich721, 1)) ./ ...
    (nLabelPermutations + 1);

PermutationTests = table(methodNames, observedAP.', apPermP.', ...
    observedEnrich721.', enrichPermP.', ...
    'VariableNames', {'Method', 'ObservedAP', 'AP_PermutationP', ...
    'ObservedEnrichmentTop721', 'EnrichmentTop721_PermutationP'});

MethodMetrics.AP_PermutationP = apPermP.';
MethodMetrics.EnrichmentTop721_PermutationP = enrichPermP.';

%% 8. Paired stratified bootstrap
fprintf('Running %d paired stratified bootstrap replicates...\n', nBootstrap);
% Stratified bootstrap samples positive and negative labels separately.
posIndex = find(isMarkT);
negIndex = find(isMarkF);
bootAP = zeros(nBootstrap, nMethod);
bootEnrich721 = zeros(nBootstrap, nMethod);

for b = 1:nBootstrap
    sampledPos = posIndex(randi(numel(posIndex), numel(posIndex), 1));
    sampledNeg = negIndex(randi(numel(negIndex), numel(negIndex), 1));
    sampledIndex = [sampledPos; sampledNeg];
    sampledLabels = [true(numel(sampledPos),1); false(numel(sampledNeg),1)];
    sampledSecondary = (1:numel(sampledIndex)).';

    for m = 1:nMethod
        sampledScores = methodScores(sampledIndex, m);
        ord = deterministicRankOrder(sampledScores, sampledSecondary);
        yRank = sampledLabels(ord);
        bootAP(b,m) = averagePrecisionFromRankedLabels(yRank);
        bootEnrich721(b,m) = mean(yRank(1:fixedTopK)) / prevalence;
    end
end

pairI = [1; 1; 2];
pairJ = [2; 3; 3];
nPair = numel(pairI);
BootstrapRows = cell(nPair * 2, 1);
rowCounter = 0;

for p = 1:nPair
    i = pairI(p);
    j = pairJ(p);

    dAP = bootAP(:,i) - bootAP(:,j);
    ciAP = empiricalCI(dAP, bootstrapAlpha);
    pAP = bootstrapTwoSidedP(dAP);
    rowCounter = rowCounter + 1;
    BootstrapRows{rowCounter} = table( ...
        methodNames(i) + " - " + methodNames(j), "Average precision", ...
        observedAP(i) - observedAP(j), ciAP(1), ciAP(2), pAP, ...
        'VariableNames', {'Comparison', 'Metric', 'ObservedDifference', ...
        'CI95_Lower', 'CI95_Upper', 'BootstrapTwoSidedP'});

    dE = bootEnrich721(:,i) - bootEnrich721(:,j);
    ciE = empiricalCI(dE, bootstrapAlpha);
    pE = bootstrapTwoSidedP(dE);
    rowCounter = rowCounter + 1;
    BootstrapRows{rowCounter} = table( ...
        methodNames(i) + " - " + methodNames(j), "Enrichment at Top 721", ...
        observedEnrich721(i) - observedEnrich721(j), ciE(1), ciE(2), pE, ...
        'VariableNames', {'Comparison', 'Metric', 'ObservedDifference', ...
        'CI95_Lower', 'CI95_Upper', 'BootstrapTwoSidedP'});
end

BootstrapComparisons = vertcat(BootstrapRows{:});

%% 9. Latent-score sensitivity analysis
SensitivityScores = struct();
SensitivityScores.HR_GGM_ExplicitLH = hrScore;

SensitivityScores.LVGGM_PriorSelectedComponent = ...
    lvScores.PriorSelectedComponent;
SensitivityScores.LVGGM_SeedContrastSubspace = lvScores.SeedContrastSubspace;
SensitivityScores.LVGGM_SeedContrast1D = lvScores.SeedContrast1D;
SensitivityScores.LVGGM_SeedPCA1D = lvScores.SeedPCA1D;
SensitivityScores.LVGGM_SeedPCA90Subspace = lvScores.SeedPCA90Subspace;
SensitivityScores.LVGGM_DominantEigenvector = lvScores.DominantEigenvector;
SensitivityScores.LVGGM_OverallLatentMagnitude = lvScores.OverallLatentMagnitude;
if runOracleDiagnostic
    SensitivityScores.LVGGM_OracleBestTopEigencomponent = ...
        lvScores.OracleBestTopEigencomponent;
end

SensitivityScores.Prior_Correlation = priorCorrScore;

sensitivityNames = string(fieldnames(SensitivityScores));
nSensitivity = numel(sensitivityNames);
SensitivityRows = cell(nSensitivity, 1);

for s = 1:nSensitivity
    fieldName = char(sensitivityNames(s));
    score = SensitivityScores.(fieldName);
    er = evaluateRanking(score, isMarkT, KList, fixedTopK, originalEvalIndex);
    isOracle = contains(sensitivityNames(s), "Oracle", 'IgnoreCase', true);
    SensitivityRows{s} = table(sensitivityNames(s), isOracle, er.AP, ...
        er.NormalizedAP, er.AUROC, er.PositiveRateAtFixedK, ...
        er.EnrichmentAtFixedK, er.OddsRatio, er.FisherRightP, ...
        'VariableNames', {'ScoreDefinition', 'UsesEvaluationLabels', ...
        'AveragePrecision', 'NormalizedAP', 'AUROC', ...
        'PositiveRateTop721', 'EnrichmentTop721', ...
        'OddsRatioTop721', 'FisherRightP'});
end
SensitivityMetrics = vertcat(SensitivityRows{:});

%% 10. Low-rank diagnostic tables
LowRankDiagnostics = lvDiag;
SpectrumLong = lvSpectrum;
SeedDirectionNull = lvDirectionNull;

%% 11. Save numerical outputs
writetable(MethodMetrics, fullfile(outputFolder, ...
    'Main_Method_Metrics.csv'));
writetable(TopKLong, fullfile(outputFolder, ...
    'Main_TopK_Results_Long.csv'));
writetable(RankedLong, fullfile(outputFolder, ...
    'Ranked_Genes_All_Methods.csv'));
writetable(ContingencyStats, fullfile(outputFolder, ...
    'Top721_Contingency_and_Fisher.csv'));
writetable(PermutationTests, fullfile(outputFolder, ...
    'Label_Permutation_Tests.csv'));
writetable(BootstrapComparisons, fullfile(outputFolder, ...
    'Paired_Bootstrap_Method_Differences.csv'));
writetable(SensitivityMetrics, fullfile(outputFolder, ...
    'Latent_Score_Sensitivity_Metrics.csv'));
writetable(LowRankDiagnostics, fullfile(outputFolder, ...
    'LowRank_Spectral_Diagnostics.csv'));
writetable(SpectrumLong, fullfile(outputFolder, ...
    'LowRank_Eigenvalue_Spectra.csv'));
writetable(SeedDirectionNull, fullfile(outputFolder, ...
    'Seed_Direction_Permutation_Null.csv'));

writetable(FigureMethodMetrics, fullfile(outputFolder, ...
    'Figure_Method_Metrics_HR_LV_PriorCorrelation.csv'));
writetable(FigureTopKLong, fullfile(outputFolder, ...
    'Figure_TopK_Results_HR_LV_PriorCorrelation.csv'));

% Save bootstrap distributions for independent verification.
BootstrapDistributions = table((1:nBootstrap).', ...
    bootAP(:,1), bootAP(:,2), bootAP(:,3), ...
    bootEnrich721(:,1), bootEnrich721(:,2), bootEnrich721(:,3), ...
    'VariableNames', {'BootstrapReplicate', ...
    'AP_HRGGM', 'AP_LVGGM', 'AP_PriorCorrelation', ...
    'Enrichment721_HRGGM', 'Enrichment721_LVGGM', ...
    'Enrichment721_PriorCorrelation'});
writetable(BootstrapDistributions, fullfile(outputFolder, ...
    'Bootstrap_Distributions.csv'));

%% 12. Generate the comparison figure
% Plot appearance and baseline settings.
colors = [0.0000, 0.4470, 0.7410; ...
          0.8500, 0.3250, 0.0980; ...
          0.4660, 0.6740, 0.1880];
lineStyles = {'-', '--', '-.'};
markers = {'o', 's', '^'};
baselineColor = [0.35 0.35 0.35];

% Create a 178 mm wide, three-panel figure.
figMain = figure('Color', 'w', 'Units', 'inches', ...
    'Position', [0.5, 0.5, 7.01, 2.65], 'PaperPositionMode', 'auto');
tl = tiledlayout(figMain, 1, 3, 'TileSpacing', 'compact', ...
    'Padding', 'compact');

% Panel A: Top-K positive rate.
ax1 = nexttile(tl, 1);
hold(ax1, 'on');
hA = gobjects(nFigureMethod + 1, 1);
for m = 1:nFigureMethod
    rows = FigureTopKLong.Method == figureMethodNames(m);
    hA(m) = plot(ax1, FigureTopKLong.TopK(rows), FigureTopKLong.PositiveRate(rows), ...
        'LineStyle', lineStyles{m}, 'Marker', markers{m}, ...
        'Color', colors(m,:), 'LineWidth', 1.0, 'MarkerSize', 3.0, ...
        'MarkerFaceColor', 'none', 'DisplayName', char(figureMethodNames(m)));
end
hA(end) = yline(ax1, prevalence, ':', ...
    'Color', baselineColor, 'LineWidth', 0.9, ...
    'DisplayName', 'Baseline');
xlabel(ax1, 'Top-K non-prior genes');
ylabel(ax1, 'Mark-T positive rate');
title(ax1, '(A) Top-K positive rate');
xlim(ax1, [0 nEval]);
ylim(ax1, [0 max(0.60, max(FigureTopKLong.PositiveRate) * 1.08)]);
grid(ax1, 'on'); box(ax1, 'on');
legend(ax1, hA, 'Location', 'best', 'Box', 'off', 'FontSize', 5.5);

% Panel B: Top-K enrichment.
ax2 = nexttile(tl, 2);
hold(ax2, 'on');
hB = gobjects(nFigureMethod + 1, 1);
for m = 1:nFigureMethod
    rows = FigureTopKLong.Method == figureMethodNames(m);
    hB(m) = plot(ax2, FigureTopKLong.TopK(rows), FigureTopKLong.EnrichmentFold(rows), ...
        'LineStyle', lineStyles{m}, 'Marker', markers{m}, ...
        'Color', colors(m,:), 'LineWidth', 1.0, 'MarkerSize', 3.0, ...
        'MarkerFaceColor', 'none', 'DisplayName', char(figureMethodNames(m)));
end
hB(end) = yline(ax2, 1, ':', ...
    'Color', baselineColor, 'LineWidth', 0.9, ...
    'DisplayName', 'Baseline');
xlabel(ax2, 'Top-K non-prior genes');
ylabel(ax2, 'Fold enrichment');
title(ax2, '(B) Mark-T enrichment');
xlim(ax2, [0 nEval]);
ylim(ax2, [0.8 max(1.2, max(FigureTopKLong.EnrichmentFold) * 1.08)]);
grid(ax2, 'on'); box(ax2, 'on');
legend(ax2, hB, 'Location', 'best', 'Box', 'off', 'FontSize', 5.5);

% Panel C: precision-recall curves.
ax3 = nexttile(tl, 3);
hold(ax3, 'on');
hC = gobjects(nFigureMethod + 1, 1);
for m = 1:nFigureMethod
    hC(m) = plot(ax3, FigurePR{m}.Recall, FigurePR{m}.Precision, ...
        'LineStyle', lineStyles{m}, 'Color', colors(m,:), ...
        'LineWidth', 1.2, 'DisplayName', char(figureMethodNames(m)));
end
hC(end) = yline(ax3, prevalence, ':', ...
    'Color', baselineColor, 'LineWidth', 0.9, ...
    'DisplayName', 'Baseline');
xlabel(ax3, 'Recall');
ylabel(ax3, 'Precision');
title(ax3, '(C) Precision-recall curves');
xlim(ax3, [0 1]);
ylim(ax3, [0 1]);
grid(ax3, 'on'); box(ax3, 'on');
legend(ax3, hC, 'Location', 'best', 'Box', 'off', 'FontSize', 5.5);

formatAxes([ax1, ax2, ax3]);
set([ax1.XLabel ax1.YLabel ax2.XLabel ax2.YLabel ax3.XLabel ax3.YLabel], ...
    'FontSize', 7.2, 'FontName', 'Arial');
set([ax1.Title ax2.Title ax3.Title], ...
    'FontSize', 7.5, 'FontName', 'Arial', 'FontWeight', 'bold');

exportgraphics(figMain, fullfile(outputFolder, ...
    'MethodComparison_Main_ABC_1x3.png'), 'Resolution', 1200);

%% 13. Write the analysis summary
summaryFile = fullfile(outputFolder, 'Comparison_Summary.txt');
fid = fopen(summaryFile, 'w');
if fid < 0
    error('Cannot create summary file: %s', summaryFile);
end
cleanupSummary = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, 'HR-GGM, LV-GGM, and prior-correlation comparison\n');
fprintf(fid, '=================================================\n\n');
fprintf(fid, 'Evaluation genes: %d\n', nEval);
fprintf(fid, 'Prior genes excluded from evaluation: %d\n', nPriorSeeds);
fprintf(fid, 'Mark-T: %d\n', nMarkT);
fprintf(fid, 'Mark-F: %d\n', nMarkF);
fprintf(fid, 'Baseline Mark-T prevalence: %.12f\n\n', prevalence);

fprintf(fid, 'Primary score definitions\n');
fprintf(fid, '-------------------------\n');
fprintf(fid, 'HR-GGM: absolute gene-hidden Ca precision weight.\n');
fprintf(fid, ['LV-GGM: one canonical eigencomponent selected exactly as in ' ...
    'the simulation analysis by the largest mean standardized loading ' ...
    'among the 42 prior genes; gene score = sqrt(lambda_k)*abs(u_k). ' ...
    'No Mark-T/F label is used for selection.\n']);
fprintf(fid, ['Prior correlation: mean absolute Pearson correlation with ' ...
    'the 42 prior genes across samples.\n\n']);

fprintf(fid, 'Main method metrics\n');
fprintf(fid, '-------------------\n');
for m = 1:nMethod
    fprintf(fid, '%s\n', char(methodNames(m)));
    fprintf(fid, '  AP: %.12f\n', MethodMetrics.AveragePrecision(m));
    fprintf(fid, '  Normalized AP: %.12f\n', MethodMetrics.NormalizedAP(m));
    fprintf(fid, '  AUROC: %.12f\n', MethodMetrics.AUROC(m));
    fprintf(fid, '  Top-721 positive rate: %.12f\n', ...
        MethodMetrics.PositiveRateTop721(m));
    fprintf(fid, '  Top-721 enrichment: %.12f\n', ...
        MethodMetrics.EnrichmentTop721(m));
    fprintf(fid, '  Top-721 OR: %.12f [%.12f, %.12f]\n', ...
        MethodMetrics.OddsRatioTop721(m), ...
        MethodMetrics.OR_CI95_Lower(m), ...
        MethodMetrics.OR_CI95_Upper(m));
    fprintf(fid, '  Fisher right-tail p: %.12e\n', ...
        MethodMetrics.FisherRightP(m));
    fprintf(fid, '  AP label-permutation p: %.12e\n', ...
        MethodMetrics.AP_PermutationP(m));
    fprintf(fid, '  Top-721 enrichment permutation p: %.12e\n\n', ...
        MethodMetrics.EnrichmentTop721_PermutationP(m));
end

fprintf(fid, 'Low-rank diagnostics\n');
fprintf(fid, '--------------------\n');
for r = 1:height(LowRankDiagnostics)
    fprintf(fid, '%s\n', char(LowRankDiagnostics.Method(r)));
    fprintf(fid, '  Input size: %d x %d\n', ...
        LowRankDiagnostics.InputRows(r), LowRankDiagnostics.InputColumns(r));
    fprintf(fid, '  Sign flipped before PSD extraction: %d\n', ...
        LowRankDiagnostics.SignFlipped(r));
    fprintf(fid, '  Numerical PSD rank: %d\n', ...
        LowRankDiagnostics.NumericalRankPSD(r));
    fprintf(fid, '  Rank explaining 90%% / 95%% / 99%% energy: %d / %d / %d\n', ...
        LowRankDiagnostics.Rank90(r), LowRankDiagnostics.Rank95(r), ...
        LowRankDiagnostics.Rank99(r));
    fprintf(fid, '  Prior-selected canonical component: %d\n', ...
        LowRankDiagnostics.PriorSelectedComponent(r));
    fprintf(fid, '  Prior enrichment of selected component: %.12f\n', ...
        LowRankDiagnostics.PriorSelectedEnrichment(r));
    fprintf(fid, '  Significant seed-enriched directions: %d\n', ...
        LowRankDiagnostics.SelectedSeedDirections(r));
    fprintf(fid, '  Seed-direction null threshold: %.12f\n\n', ...
        LowRankDiagnostics.SeedDirectionNullThreshold(r));
end

fprintf(fid, 'Interpretation guardrails\n');
fprintf(fid, '-------------------------\n');
fprintf(fid, ['1. Rank greater than one is not itself a failure; it indicates ' ...
    'multiple latent directions.\n']);
fprintf(fid, ['2. The primary LV-GGM score uses one prior-selected canonical ' ...
    'eigencomponent selected from the prior-gene loadings.\n']);
fprintf(fid, ['3. The independent Mark-T/F labels are used only for final ' ...
    'evaluation, permutation tests, bootstrap comparisons, and the explicitly ' ...
    'labeled oracle diagnostic.\n']);
fprintf(fid, ['4. Oracle rows in Latent_Score_Sensitivity_Metrics.csv are ' ...
    'descriptive upper bounds and must not be reported as inferential results.\n']);

fprintf('\n============================================================\n');
fprintf('Analysis completed successfully.\n');
fprintf('All outputs were written to: %s\n', outputFolder);
fprintf('============================================================\n');
disp(MethodMetrics);
disp(FigureMethodMetrics);
disp(LowRankDiagnostics);
disp(BootstrapComparisons);

%% Local functions
function filePath = resolveExistingFile(candidateFiles)
% Return the first existing path from a list of candidate filenames.
    for i = 1:numel(candidateFiles)
        if exist(candidateFiles{i}, 'file') == 2
            filePath = candidateFiles{i};
            return;
        end
    end
    error('None of the required files exists: %s', strjoin(candidateFiles, ', '));
end

function variableName = findVariableName(T, candidateNames)
% Match a table variable name after normalizing header text.
    available = string(T.Properties.VariableNames);
    normalizedAvailable = normalizeHeader(available);
    normalizedCandidates = normalizeHeader(string(candidateNames));
    variableName = '';
    for i = 1:numel(normalizedCandidates)
        idx = find(normalizedAvailable == normalizedCandidates(i), 1);
        if ~isempty(idx)
            variableName = T.Properties.VariableNames{idx};
            return;
        end
    end
    error('Required column not found. Available columns: %s', ...
        strjoin(T.Properties.VariableNames, ', '));
end

function normalized = normalizeHeader(x)
% Normalize table headers for case- and punctuation-insensitive matching.
    normalized = lower(regexprep(string(x), '[^A-Za-z0-9]', ''));
end

function x = numericColumn(rawColumn)
% Convert a table column to a numeric column vector.
    if isnumeric(rawColumn)
        x = double(rawColumn);
    else
        x = str2double(string(rawColumn));
    end
    x = x(:);
end

function [scores, diagTable, spectrumTable, nullTable] = ...
    processLatentMatrix(filePath, methodName, opt)
% Validate a latent matrix, extract spectral scores, and return diagnostics.

    fprintf('  Reading %s...\n', filePath);
    Lraw = readmatrix(filePath);
    inputRows = size(Lraw,1);
    inputCols = size(Lraw,2);

    if inputRows ~= inputCols
        error('%s latent matrix is not square: %d x %d.', ...
            methodName, inputRows, inputCols);
    end
    if inputRows == opt.nGene
        L = Lraw;
    elseif inputRows == opt.nGene + 2
        fprintf('  %s is 3002 x 3002; using the first 3000 gene rows/columns.\n', ...
            methodName);
        L = Lraw(1:opt.nGene, 1:opt.nGene);
    else
        error(['%s latent matrix must be 3000 x 3000 or 3002 x 3002; ' ...
               'found %d x %d.'], methodName, inputRows, inputCols);
    end
    clear Lraw;

    if any(~isfinite(L), 'all')
        error('%s latent matrix contains NaN or Inf.', methodName);
    end

    symmetryError = max(abs(L - L.'), [], 'all');
    L = (L + L.') / 2;
    froNorm = norm(L, 'fro');
    fprintf('  Symmetry error before symmetrization: %.3e\n', symmetryError);
    fprintf('  Computing full symmetric eigendecomposition...\n');

    [U, lambda] = eig(L, 'vector');
    lambda = real(lambda);
    U = real(U);
    [lambda, ord] = sort(lambda, 'descend');
    U = U(:,ord);
    clear L;

    positiveEnergyOriginal = sum(max(lambda, 0));
    negativeEnergyOriginal = sum(max(-lambda, 0));
    signFlipped = negativeEnergyOriginal > positiveEnergyOriginal;
    if signFlipped
        fprintf(['  Negative spectral energy dominates; multiplying the ' ...
                 'latent matrix by -1 before PSD extraction.\n']);
        lambda = -lambda;
        [lambda, ord2] = sort(lambda, 'descend');
        U = U(:,ord2);
    end

    maxLambda = max(lambda);
    eigTol = max(opt.absoluteEigenTolerance, ...
        opt.relativeEigenTolerance * max(maxLambda, eps));
    keep = lambda > eigTol;
    lambdaPos = lambda(keep);
    Upos = U(:,keep);
    clear U;

    numericalRank = numel(lambdaPos);
    if numericalRank == 0
        error('%s has no positive eigenvalue above tolerance %.3e.', ...
            methodName, eigTol);
    end

    positiveEnergy = sum(lambdaPos);
    residualNegativeEnergy = sum(max(-lambda, 0));
    totalAbsEnergy = sum(abs(lambda));
    residualNegativeFraction = residualNegativeEnergy / max(totalAbsEnergy, eps);

    cumulativeEnergy = cumsum(lambdaPos) / positiveEnergy;
    rank90 = find(cumulativeEnergy >= 0.90, 1, 'first');
    rank95 = find(cumulativeEnergy >= 0.95, 1, 'first');
    rank99 = find(cumulativeEnergy >= 0.99, 1, 'first');
    top1Share = lambdaPos(1) / positiveEnergy;

    fprintf('  Numerical PSD rank: %d\n', numericalRank);
    fprintf('  Top eigenvalue energy share: %.6f\n', top1Share);
    fprintf('  Rank for 90%% / 95%% / 99%% energy: %d / %d / %d\n', ...
        rank90, rank95, rank99);

    B = Upos .* sqrt(lambdaPos.');
    clear Upos;

    canonicalFeatureScores = abs(B);
    priorEnrichment = nan(size(canonicalFeatureScores, 2), 1);

    for k = 1:size(canonicalFeatureScores, 2)
        scoreZ = zscoreSafeVector(canonicalFeatureScores(:, k));
        priorEnrichment(k) = mean(scoreZ(opt.seedIndex), 'omitnan');
    end

    priorEnrichment(~isfinite(priorEnrichment)) = -Inf;
    [priorSelectedEnrichment, priorSelectedComponent] = ...
        max(priorEnrichment);

    if ~isfinite(priorSelectedEnrichment)
        error('%s has no finite prior-enrichment score.', methodName);
    end

    scorePriorSelectedAll = ...
        canonicalFeatureScores(:, priorSelectedComponent);

    fprintf('  Prior-selected canonical component: %d\n', ...
        priorSelectedComponent);
    fprintf('  Mean standardized prior loading: %.6f\n', ...
        priorSelectedEnrichment);

    seedRows = opt.seedIndex;
    evalRows = opt.evalIndex;
    Bseed = B(seedRows,:);
    nSeed = numel(seedRows);
    nAll = size(B,1);
    r = size(B,2);

    [~, Sseed, VseedContrast] = svd(Bseed, 'econ');
    contrastEigenvalues = (diag(Sseed).^2) / nSeed;

    nullMax = zeros(opt.nSeedPermutations, 1);
    for b = 1:opt.nSeedPermutations
        randomSeedRows = randperm(nAll, nSeed);
        Bperm = B(randomSeedRows,:);
        Gram = (Bperm * Bperm.') / nSeed;
        Gram = (Gram + Gram.') / 2;
        nullMax(b) = max(real(eig(Gram)));
    end
    nullThreshold = empiricalQuantile(nullMax, 1 - opt.seedDirectionAlpha);

    selected = find(contrastEigenvalues > nullThreshold);
    noSignificantDirection = isempty(selected);
    if noSignificantDirection
        selected = 1;
        fprintf(['  No raw seed-enriched direction exceeded the %.1f%% FWER ' ...
                 'threshold; the top raw seed direction is retained only ' ...
                 'to permit a conservative ranking comparison.\n'], ...
                 100 * (1 - opt.seedDirectionAlpha));
    end
    selected = selected(1:min(numel(selected), opt.maxSeedDirections));
    nSelected = numel(selected);

    Qcontrast = VseedContrast(:,selected);
    projectionContrast = B * Qcontrast;
    scoreContrastSubspaceAll = sqrt(sum(projectionContrast.^2, 2));
    scoreContrast1DAll = abs(projectionContrast(:,1));

    fprintf('  Latent-dimension whitening: OFF\n');
    fprintf('  Raw seed-direction null threshold: %.6f\n', nullThreshold);
    fprintf('  Retained raw seed-enriched directions: %d\n', nSelected);

    [~, SseedRaw, VseedRaw] = svd(Bseed, 'econ');
    scoreSeedPCA1DAll = abs(B * VseedRaw(:,1));

    seedSingularEnergy = diag(SseedRaw).^2;
    if sum(seedSingularEnergy) > 0
        cumulativeSeedEnergy = cumsum(seedSingularEnergy) / ...
            sum(seedSingularEnergy);
        seedPcaDim = find(cumulativeSeedEnergy >= opt.seedPcaEnergyTarget, ...
            1, 'first');
    else
        seedPcaDim = 1;
    end
    seedPcaDim = max(seedPcaDim, 1);
    scoreSeedPCA90All = sqrt(sum((B * VseedRaw(:,1:seedPcaDim)).^2, 2));

    scoreDominantAll = abs(B(:,1));
    scoreMagnitudeAll = sqrt(sum(B.^2, 2));

    oracleComponent = NaN;
    oracleAP = NaN;
    scoreOracleAll = nan(nAll,1);
    if opt.runOracleDiagnostic
        nOracle = min([opt.oracleMaxComponents, r, rank99]);
        bestAP = -inf;
        bestJ = 1;
        for j = 1:nOracle
            candidate = abs(B(evalRows,j));
            ordCandidate = deterministicRankOrder(candidate, ...
                opt.originalEvalIndex);
            apCandidate = averagePrecisionFromRankedLabels( ...
                opt.evalLabels(ordCandidate));
            if apCandidate > bestAP
                bestAP = apCandidate;
                bestJ = j;
            end
        end
        oracleComponent = bestJ;
        oracleAP = bestAP;
        scoreOracleAll = abs(B(:,bestJ));
    end

    scores = struct();
    scores.PriorSelectedComponent = scorePriorSelectedAll(evalRows);
    scores.SeedContrastSubspace = scoreContrastSubspaceAll(evalRows);
    scores.SeedContrast1D = scoreContrast1DAll(evalRows);
    scores.SeedPCA1D = scoreSeedPCA1DAll(evalRows);
    scores.SeedPCA90Subspace = scoreSeedPCA90All(evalRows);
    scores.DominantEigenvector = scoreDominantAll(evalRows);
    scores.OverallLatentMagnitude = scoreMagnitudeAll(evalRows);
    scores.OracleBestTopEigencomponent = scoreOracleAll(evalRows);

    diagTable = table(string(methodName), inputRows, inputCols, ...
        symmetryError, froNorm, signFlipped, positiveEnergyOriginal, ...
        negativeEnergyOriginal, residualNegativeFraction, eigTol, ...
        numericalRank, top1Share, rank90, rank95, rank99, ...
        priorSelectedComponent, priorSelectedEnrichment, ...
        seedPcaDim, nullThreshold, nSelected, noSignificantDirection, ...
        oracleComponent, oracleAP, ...
        'VariableNames', {'Method', 'InputRows', 'InputColumns', ...
        'SymmetryError', 'FrobeniusNorm', 'SignFlipped', ...
        'PositiveEnergyBeforeOrientation', ...
        'NegativeEnergyBeforeOrientation', ...
        'ResidualNegativeEnergyFraction', 'EigenTolerance', ...
        'NumericalRankPSD', 'Top1EnergyShare', 'Rank90', 'Rank95', ...
        'Rank99', 'PriorSelectedComponent', 'PriorSelectedEnrichment', ...
        'SeedPCA90Dimension', 'SeedDirectionNullThreshold', ...
        'SelectedSeedDirections', 'NoSignificantSeedDirection', ...
        'OracleBestComponent', 'OracleBestComponentAP'});

    componentRank = (1:numericalRank).';
    spectrumTable = table(repmat(string(methodName), numericalRank, 1), ...
        componentRank, lambdaPos, lambdaPos / positiveEnergy, ...
        cumulativeEnergy, ...
        'VariableNames', {'Method', 'ComponentRank', 'Eigenvalue', ...
        'EnergyFraction', 'CumulativeEnergy'});

    nullTable = table(repmat(string(methodName), opt.nSeedPermutations, 1), ...
        (1:opt.nSeedPermutations).', nullMax, ...
        repmat(nullThreshold, opt.nSeedPermutations, 1), ...
        'VariableNames', {'Method', 'Permutation', ...
        'MaximumSeedContrastEigenvalue', 'FWERThreshold'});
end

function result = evaluateRanking(score, y, KList, fixedK, secondaryIndex)
% Compute ranking, Top-K, PR, ROC, and contingency-based metrics.
    score = double(score(:));
    y = logical(y(:));
    secondaryIndex = secondaryIndex(:);
    n = numel(y);
    nPos = sum(y);
    prevalence = nPos / n;

    if numel(score) ~= n || numel(secondaryIndex) ~= n
        error('Score, labels, and secondary index must have equal length.');
    end
    if any(~isfinite(score))
        error('Ranking score contains NaN or Inf.');
    end
    if any(KList < 1 | KList > n) || fixedK < 1 || fixedK > n
        error('Invalid K value.');
    end

    ord = deterministicRankOrder(score, secondaryIndex);
    yRank = y(ord);
    scoreRank = score(ord);
    cumulativePos = cumsum(yRank);
    rankVector = (1:n).';

    runningPositiveRate = cumulativePos ./ rankVector;
    runningRecall = cumulativePos / nPos;
    runningEnrichment = runningPositiveRate / prevalence;

    positiveCount = cumulativePos(KList);
    positiveRate = positiveCount ./ KList;
    recall = positiveCount / nPos;
    enrichment = positiveRate / prevalence;

    AP = averagePrecisionFromRankedLabels(yRank);
    normalizedAP = AP / prevalence;
    [rocFpr, rocTpr, auc] = rocCurveTieAware(score, y);
    [prRecall, prPrecision] = prCurveFromRankedLabels(yRank);

    modelPositive = false(n,1);
    modelPositive(ord(1:fixedK)) = true;
    a = sum(y & modelPositive);
    b = sum(y & ~modelPositive);
    c = sum(~y & modelPositive);
    d = sum(~y & ~modelPositive);
    contingency = [a b; c d];

    [OR, ci95] = oddsRatioWithCI(a,b,c,d);
    fisherP = fisherExactRightTail(a,b,c,d);

    cutoffScore = scoreRank(fixedK);
    if fixedK < n
        nextScore = scoreRank(fixedK + 1);
    else
        nextScore = NaN;
    end
    tieMask = scoreRank == cutoffScore;
    tiedAtCutoff = sum(tieMask);

    result = struct();
    result.RankOrder = ord;
    result.RunningPositiveRate = runningPositiveRate;
    result.RunningRecall = runningRecall;
    result.RunningEnrichment = runningEnrichment;
    result.TopKPositiveCount = positiveCount;
    result.TopKPositiveRate = positiveRate;
    result.TopKRecall = recall;
    result.TopKEnrichment = enrichment;
    result.AP = AP;
    result.NormalizedAP = normalizedAP;
    result.AUROC = auc;
    result.PositiveRateAtFixedK = positiveRate(KList == fixedK);
    if isempty(result.PositiveRateAtFixedK)
        result.PositiveRateAtFixedK = cumulativePos(fixedK) / fixedK;
    end
    result.EnrichmentAtFixedK = result.PositiveRateAtFixedK / prevalence;
    result.RecallAtFixedK = cumulativePos(fixedK) / nPos;
    result.Contingency = contingency;
    result.OddsRatio = OR;
    result.OR_CI95 = ci95;
    result.FisherRightP = fisherP;
    result.CutoffScore = cutoffScore;
    result.NextScore = nextScore;
    result.TiedAtCutoff = tiedAtCutoff;
    result.PR = table(prRecall, prPrecision, ...
        'VariableNames', {'Recall', 'Precision'});
    result.ROC = table(rocFpr, rocTpr, ...
        'VariableNames', {'FPR', 'TPR'});
end

function ord = deterministicRankOrder(score, secondaryIndex)
% Sort scores in descending order and resolve ties by the original index.
    score = score(:);
    secondaryIndex = secondaryIndex(:);
    [~, ord] = sortrows([-score, secondaryIndex], [1 2]);
end

function ap = averagePrecisionFromRankedLabels(yRank)
% Compute non-interpolated average precision from ranked binary labels.
    yRank = logical(yRank(:));
    nPos = sum(yRank);
    if nPos == 0
        ap = NaN;
        return;
    end
    precisionAtRank = cumsum(yRank) ./ (1:numel(yRank)).';
    ap = sum(precisionAtRank(yRank)) / nPos;
end

function [recall, precision] = prCurveFromRankedLabels(yRank)
% Generate precision-recall coordinates from ranked binary labels.
    yRank = logical(yRank(:));
    nPos = sum(yRank);
    tp = cumsum(yRank);
    precision = tp ./ (1:numel(yRank)).';
    recall = tp / nPos;
    recall = [0; recall];
    precision = [1; precision];
end

function [fpr, tpr, auc] = rocCurveTieAware(score, y)
% Generate a tie-aware ROC curve and calculate its area.
    score = double(score(:));
    y = logical(y(:));
    [scoreSorted, ord] = sort(score, 'descend');
    ySorted = y(ord);

    groupEnd = [find(diff(scoreSorted) ~= 0); numel(scoreSorted)];
    tpCum = cumsum(ySorted);
    fpCum = cumsum(~ySorted);
    nPos = sum(y);
    nNeg = sum(~y);

    tpr = [0; tpCum(groupEnd) / nPos];
    fpr = [0; fpCum(groupEnd) / nNeg];
    auc = trapz(fpr, tpr);
end

function [OR, ci95] = oddsRatioWithCI(a,b,c,d)
% Calculate an odds ratio and Wald confidence interval.
    cells = double([a b c d]);
    if any(cells == 0)
        cells = cells + 0.5;
    end
    aa = cells(1); bb = cells(2); cc = cells(3); dd = cells(4);
    OR = (aa * dd) / (bb * cc);
    se = sqrt(1/aa + 1/bb + 1/cc + 1/dd);
    ci95 = exp(log(OR) + [-1 1] * 1.96 * se);
end

function p = fisherExactRightTail(a,b,c,d)
% Calculate the right-tailed Fisher exact-test probability.
    rowPos = a + b;
    rowNeg = c + d;
    colSelected = a + c;
    totalN = a + b + c + d;

    maxX = min(rowPos, colSelected);
    x = a:maxX;
    logP = logChoose(rowPos, x) + ...
           logChoose(rowNeg, colSelected - x) - ...
           logChoose(totalN, colSelected);
    maxLogP = max(logP);
    p = exp(maxLogP) * sum(exp(logP - maxLogP));
    p = min(max(p, 0), 1);
end

function y = logChoose(n,k)
% Calculate log binomial coefficients for vector inputs.
    y = -inf(size(k));
    valid = k >= 0 & k <= n & floor(k) == k;
    kv = k(valid);
    y(valid) = gammaln(n + 1) - gammaln(kv + 1) - ...
        gammaln(n - kv + 1);
end

function ci = empiricalCI(x, alpha)
% Calculate an equal-tailed empirical confidence interval.
    ci = [empiricalQuantile(x, alpha/2), ...
          empiricalQuantile(x, 1 - alpha/2)];
end

function q = empiricalQuantile(x, probability)
% Calculate a linearly interpolated empirical quantile.
    x = sort(double(x(:)));
    n = numel(x);
    if n == 0
        q = NaN;
        return;
    end
    probability = min(max(probability, 0), 1);
    position = 1 + (n - 1) * probability;
    lo = floor(position);
    hi = ceil(position);
    if lo == hi
        q = x(lo);
    else
        q = x(lo) + (position - lo) * (x(hi) - x(lo));
    end
end

function p = bootstrapTwoSidedP(differenceSamples)
% Calculate a two-sided sign-based bootstrap probability.
    x = differenceSamples(:);
    n = numel(x);
    pLower = (1 + sum(x <= 0)) / (n + 1);
    pUpper = (1 + sum(x >= 0)) / (n + 1);
    p = min(1, 2 * min(pLower, pUpper));
end

function priorScore = computePriorCorrelationScore(filePath, nGene, nPriorSeeds)
% Score each non-prior gene by mean absolute correlation with prior genes.

    X = readmatrix(filePath);
    if any(~isfinite(X), 'all')
        error('Expression matrix for Prior correlation contains NaN or Inf.');
    end

    [nRow, nCol] = size(X);
    if nRow == nGene
        geneByCell = X;
    elseif nCol == nGene
        geneByCell = X.';
    elseif nRow == nGene + 2
        geneByCell = X(1:nGene, :);
    elseif nCol == nGene + 2
        geneByCell = X(:, 1:nGene).';
    else
        error(['Expression matrix for Prior correlation must be %d x N, N x %d, ' ...
               '%d x N, or N x %d. Found %d x %d.'], ...
              nGene, nGene, nGene + 2, nGene + 2, nRow, nCol);
    end

    seedExpr = geneByCell(1:nPriorSeeds, :).';
    evalExpr = geneByCell(nPriorSeeds + 1:nGene, :).';

    C = corr(evalExpr, seedExpr, 'Rows', 'pairwise');
    priorScore = mean(abs(C), 2, 'omitnan');
    priorScore = priorScore(:);

    if any(~isfinite(priorScore))
        error('Prior-correlation scores contain NaN or Inf. Check the expression matrix.');
    end

    fprintf('  Prior-correlation baseline computed from %d genes and %d samples/cells.\n', ...
        size(geneByCell,1), size(geneByCell,2));
end

function z = zscoreSafeVector(x)
% Standardize a vector while handling zero or invalid variance.
    x = double(x(:));
    mu = mean(x, 'omitnan');
    sigma = std(x, 0, 'omitnan');

    if ~isfinite(sigma) || sigma < eps
        z = zeros(size(x));
    else
        z = (x - mu) / sigma;
    end
end

function formatAxes(axArray)
% Apply consistent publication-format axis properties.
    for i = 1:numel(axArray)
        axArray(i).FontName = 'Arial';
        axArray(i).FontSize = 6.5;
        axArray(i).LineWidth = 0.8;
        axArray(i).TickDir = 'out';
        axArray(i).GridAlpha = 0.18;
        axArray(i).MinorGridAlpha = 0.10;
        axArray(i).Layer = 'top';
    end
end
