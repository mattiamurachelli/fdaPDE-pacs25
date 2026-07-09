%% ============================================================
%% SCRIPT TO PLOT N-SIMULATION COMPARISON
%% OUR ALM VS NLOPT AUGLAG
%% ============================================================

clear;
close all;
clc;

%% ============================================================
%% CREATE PLOTS FOLDER
%% ============================================================

folder = '../Plots';

if ~exist(folder, 'dir')
    mkdir(folder);
end

%% ============================================================
%% IMPORT CSV FILES
%% ============================================================

ntest = 1;


data_OUR = readtable( ...
    strcat('../sqp_test', num2str(ntest), '.csv'), ...
    'ReadVariableNames', true, ...
    'Delimiter', ',');

data_NLOPT = readtable( ...
    strcat('../sqp_test', num2str(ntest), '_nlopt.csv'), ...
    'ReadVariableNames', true, ...
    'Delimiter', ',');

%% ============================================================
%% NORMALIZE TEXT COLUMNS
%% ============================================================

data_OUR.OptId = string(data_OUR.OptId);
data_NLOPT.OptId = string(data_NLOPT.OptId);

if ismember('Status', data_OUR.Properties.VariableNames)
    data_OUR.Status = string(data_OUR.Status);
end

if ismember('Status', data_NLOPT.Properties.VariableNames)
    data_NLOPT.Status = string(data_NLOPT.Status);
end

%% ============================================================
%% SORT BY TEST ID
%% ============================================================

data_OUR = sortrows(data_OUR, 'TestId');
data_NLOPT = sortrows(data_NLOPT, 'TestId');

%% ============================================================
%% CONVERGENCE FILTER
%% ============================================================
% Same idea as in the report: remove failed simulations before
% comparing quantitative metrics.

tol_dist = 1e-3;
tol_res  = 1e-6;

OUR_converged = ...
    data_OUR.DistSol <= tol_dist & data_OUR.Residual <= tol_res;

NLOPT_converged = ...
    data_NLOPT.DistSol <= tol_dist & data_NLOPT.Residual <= tol_res;

OUR_failures = sum(~OUR_converged);
NLOPT_failures = sum(~NLOPT_converged);

data_OUR_conv = data_OUR(OUR_converged, :);
data_NLOPT_conv = data_NLOPT(NLOPT_converged, :);

%% ============================================================
%% PLOT 1 : NUMBER OF FAILED SIMULATIONS
%% ============================================================

FailurePlot = figure();

failure_counts = [
    OUR_failures;
    NLOPT_failures
];

bar(failure_counts);

set(gca, ...
    'XTickLabel', {'OUR-ALM', 'NLOPT-AUGLAG'});

ylabel('Number of failed simulations');

title('Non-converged simulations');

grid on;

saveas(FailurePlot, ...
       fullfile(folder, strcat('SQPTest', num2str(ntest), '_FailureHistogram.png')));

%% ============================================================
%% PLOT 2 : NUMBER OF ITERATIONS / EVALUATIONS
%% ============================================================
% For OUR_ALM this is the number of ALM outer subproblems.
% For NLOPT_AUGLAG this is the number of objective evaluations
% returned by get_numevals().

IterPlot = figure();

Iter_all = [
    data_OUR_conv.NumIterOrEval;
    data_NLOPT_conv.NumIterOrEval
];

group_iter = [
    repmat({'OUR-SQP'}, height(data_OUR_conv), 1);
    repmat({'NLOPT-SQP'}, height(data_NLOPT_conv), 1)
];

boxplot(Iter_all, group_iter);

ylabel('Iterations / function evaluations');

title('Comparison of iterations/evaluations');

grid on;

saveas(IterPlot, ...
       fullfile(folder, strcat('SQPTest', num2str(ntest), '_IterEvalComparison.png')));

%% ============================================================
%% PLOT 3 : COMPUTATION TIME
%% ============================================================
% Convert from nanoseconds to milliseconds.

CompTimePlot = figure();

OUR_CompTime = data_OUR_conv.CompTime / 1e6;
NLOPT_CompTime = data_NLOPT_conv.CompTime / 1e6;

CompTime_all = [
    OUR_CompTime;
    NLOPT_CompTime
];

group_time = [
    repmat({'OUR-SQP'}, height(data_OUR_conv), 1);
    repmat({'NLOPT-SQP'}, height(data_NLOPT_conv), 1)
];

boxplot(CompTime_all, group_time);

ylabel('Computation time [ms]');

title('Comparison of computation time');

grid on;

saveas(CompTimePlot, ...
       fullfile(folder, strcat('SQPTest', num2str(ntest), '_CompTimeComparison.png')));

%% ============================================================
%% PLOT 4 : DISTANCE BETWEEN SOLUTIONS
%% ONLY KEEP TESTS WHERE BOTH METHODS CONVERGED
%% ============================================================

DistBetweenPlot = figure();

common_test_ids = intersect(data_OUR.TestId, data_NLOPT.TestId);

DistBetween = [];

for k = 1:length(common_test_ids)

    id = common_test_ids(k);

    row_OUR = find(data_OUR.TestId == id, 1);
    row_NLOPT = find(data_NLOPT.TestId == id, 1);

    if isempty(row_OUR) || isempty(row_NLOPT)
        continue;
    end

    our_ok = ...
        data_OUR.DistSol(row_OUR) <= tol_dist && ...
        data_OUR.Residual(row_OUR) <= tol_res;

    nlopt_ok = ...
        data_NLOPT.DistSol(row_NLOPT) <= tol_dist && ...
        data_NLOPT.Residual(row_NLOPT) <= tol_res;

    if our_ok && nlopt_ok

        p_OUR = parsePointLocal(data_OUR.MinPoint(row_OUR));
        p_NLOPT = parsePointLocal(data_NLOPT.MinPoint(row_NLOPT));

        DistBetween(end + 1, 1) = norm(p_OUR - p_NLOPT);

    end
end

group_between = ...
    repmat({'OUR-SQP vs NLOPT-SQP'}, length(DistBetween), 1);

boxplot(DistBetween, group_between);

ylabel('Distance between numerical solutions');

title('Comparison of distance between solutions');

grid on;

saveas(DistBetweenPlot, ...
       fullfile(folder, strcat('SQPTest', num2str(ntest), '_DistBetweenComparison.png')));

%% ============================================================
%% PRINT SUMMARY
%% ============================================================

fprintf('\n============================================================\n');
fprintf('OUR ALM vs NLOPT AUGLAG comparison\n');
fprintf('============================================================\n');
fprintf('Number of OUR_ALM simulations          : %d\n', height(data_OUR));
fprintf('Number of NLOPT_AUGLAG simulations    : %d\n', height(data_NLOPT));
fprintf('OUR_ALM failures                      : %d\n', OUR_failures);
fprintf('NLOPT_AUGLAG failures                 : %d\n', NLOPT_failures);
fprintf('Tests where both methods converged    : %d\n', length(DistBetween));
fprintf('Mean distance between final solutions : %.4e\n', mean(DistBetween));
fprintf('Max distance between final solutions  : %.4e\n', max(DistBetween));
fprintf('============================================================\n\n');

%% ============================================================
%% LOCAL FUNCTION
%% ============================================================

function p = parsePointLocal(pointString)

    pointString = string(pointString);

    parts = split(pointString, ';');

    p = str2double(parts).';

end