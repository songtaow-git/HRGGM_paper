clear;
clc;

%% Analysis settings
% All relative paths are resolved from the directory containing this script.
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end

% Input workbook and worksheet containing independent Mark-T/F annotations
% and threshold-based Model-T/F labels.
filename = fullfile(script_dir, ...
    'S5_Ca_associated_reference_label.xlsx');
sheet_name = 'Annotation';

% Columns containing the reference annotations and model labels.
GT_column = 'Mark T/F';
Pred_column = 'Model T/F';

% Binary T/F labels used in both columns.
positive_label = "T";
negative_label = "F";

%% Read and validate labels
T = readtable(filename, ...
    'Sheet', sheet_name, ...
    'VariableNamingRule', 'preserve');
T.("Model T/F") = T.("Mark T/F");
T.("Model T/F")(T.AbsWeight > 0.038) = cellstr("T");
T.("Model T/F")(T.AbsWeight <= 0.038) = cellstr("F");
required_columns = {GT_column, Pred_column};
if ~all(ismember(required_columns, T.Properties.VariableNames))
    error('The input table does not contain the required label columns.');
end

GT = strtrim(upper(string(T.(GT_column))));
Pred = strtrim(upper(string(T.(Pred_column))));

valid_GT = ismember(GT, [positive_label, negative_label]);
valid_Pred = ismember(Pred, [positive_label, negative_label]);
valid_rows = valid_GT & valid_Pred;

GT = GT(valid_rows);
Pred = Pred(valid_rows);

if isempty(GT)
    error('No rows contain valid binary labels in both columns.');
end

%% Two-by-two contingency table
TP = sum(GT == positive_label & Pred == positive_label);
FP = sum(GT == negative_label & Pred == positive_label);
FN = sum(GT == positive_label & Pred == negative_label);
TN = sum(GT == negative_label & Pred == negative_label);

total_n = TP + FP + FN + TN;
contingency_table = [TP, FN; FP, TN];

%% Descriptive classification metrics
precision = TP / (TP + FP);
recall = TP / (TP + FN);
specificity = TN / (TN + FP);
accuracy = (TP + TN) / total_n;
F1 = 2 * precision * recall / (precision + recall);

%% Enrichment statistics
% The baseline positive rate is the overall Mark-T prevalence. The Model-T
% positive rate is the Mark-T fraction among Model-T genes. Their ratio is
% reported as fold enrichment.
baseline_positive_rate = (TP + FN) / total_n;
modelT_positive_rate = TP / (TP + FP);
enrichment_fold = modelT_positive_rate / baseline_positive_rate;

%% Odds ratio and confidence interval
% A Haldane-Anscombe correction is applied only when at least one cell of
% the 2-by-2 table is zero.
if any([TP, FP, FN, TN] == 0)
    TP2 = TP + 0.5;
    FP2 = FP + 0.5;
    FN2 = FN + 0.5;
    TN2 = TN + 0.5;
else
    TP2 = TP;
    FP2 = FP;
    FN2 = FN;
    TN2 = TN;
end

odds_ratio = (TP2 * TN2) / (FP2 * FN2);
log_OR = log(odds_ratio);
SE_log_OR = sqrt(1 / TP2 + 1 / FP2 + 1 / FN2 + 1 / TN2);
CI_low = exp(log_OR - 1.96 * SE_log_OR);
CI_high = exp(log_OR + 1.96 * SE_log_OR);

%% Fisher's exact test
% Rows represent Mark-T and Mark-F genes; columns represent Model-T and
% Model-F genes.
p_fisher = NaN;
fisher_odds_ratio = NaN;

if exist('fishertest', 'file') == 2
    [~, p_fisher, stats_fisher] = fishertest(contingency_table);
    fisher_odds_ratio = stats_fisher.OddsRatio;
else
    warning(['fishertest is unavailable. Fisher''s exact test requires ', ...
        'Statistics and Machine Learning Toolbox.']);
end


%% Key statistical output
fprintf('\n=== Odds ratio and 95%% confidence interval ===\n');
fprintf('Odds ratio: %.4f\n', odds_ratio);
fprintf('95%% CI: %.4f - %.4f\n', CI_low, CI_high);

if isfinite(p_fisher)
    fprintf('\n=== Fisher exact test ===\n');
    fprintf('Fisher exact test p-value: %.4e\n', p_fisher);
    fprintf('Fisher odds ratio: %.4f\n', fisher_odds_ratio);
end

%% Chi-square test of independence
row_sums = sum(contingency_table, 2);
col_sums = sum(contingency_table, 1);
expected = row_sums * col_sums / total_n;

chi2_stat = sum((contingency_table - expected).^2 ./ expected, 'all');
df = 1;
p_chi2 = 1 - chi2cdf(chi2_stat, df);

%% Results
contingency_table_display = array2table(contingency_table, ...
    'VariableNames', {'Model_Positive', 'Model_Negative'}, ...
    'RowNames', {'Reference_Positive', 'Reference_Negative'});

summary_table = table( ...
    TP, FP, FN, TN, ...
    precision, recall, F1, accuracy, specificity, ...
    baseline_positive_rate, modelT_positive_rate, enrichment_fold, ...
    odds_ratio, CI_low, CI_high, ...
    p_fisher, fisher_odds_ratio, chi2_stat, p_chi2, ...
    'VariableNames', { ...
        'TP', 'FP', 'FN', 'TN', ...
        'Precision', 'Recall', 'F1', 'Accuracy', 'Specificity', ...
        'Baseline_Positive_Rate', 'ModelT_Positive_Rate', 'Enrichment_Fold', ...
        'Odds_Ratio', 'OR_95CI_Low', 'OR_95CI_High', ...
        'Fisher_p', 'Fisher_Odds_Ratio', 'ChiSquare_Statistic', 'ChiSquare_p' ...
    });

disp(contingency_table_display);
disp(summary_table);
