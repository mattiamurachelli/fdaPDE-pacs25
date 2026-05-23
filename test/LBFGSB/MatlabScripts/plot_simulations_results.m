%% ============================================================
%% SCRIPT TO PLOT THE RESULTS OF THE SIMULATIONS
%% OUR LBFGSB VS MIT LBFGSB OPTIMIZER
%% ============================================================

%% ============================================================
%% CREATE PLOTS FOLDER
%% ============================================================

folder = '../Plots';

if ~exist(folder, 'dir')
    mkdir(folder)
end

%% ============================================================
%% IMPORT CSV FILE
%% ============================================================

data = readtable( ...
    '../simulations_rosenbrock.csv', ...
    'ReadVariableNames', true, ...
    'Delimiter', ',');

%% ============================================================
%% SEPARATE METHODS
%% ============================================================

data_MIT = data(strcmp(data.OptId, 'MIT'), :);

data_OUR = data(strcmp(data.OptId, 'OUR'), :);

%% ============================================================
%% REMOVE NON-CONVERGED SIMULATIONS
%% ============================================================

tol = 1e-3;

MIT_failures = sum(data_MIT.DistSol > tol);

OUR_failures = sum(data_OUR.DistSol > tol);

% Keep only converged simulations
data_MIT = data_MIT(data_MIT.DistSol <= tol, :);

data_OUR = data_OUR(data_OUR.DistSol <= tol, :);

%% ============================================================
%% IMPORT DATA
%% ============================================================

%% Number of iterations

MIT_NumIter = data_MIT.NumIter;

OUR_NumIter = data_OUR.NumIter;

%% Computation time
%% Convert from ns to ms

MIT_CompTime = data_MIT.CompTime / (1e6);

OUR_CompTime = data_OUR.CompTime / (1e6);

%% Minimum points

MIT_MinPoint = parsePoints(data_MIT.MinPoint);

OUR_MinPoint = parsePoints(data_OUR.MinPoint);

%% Objective function value

MIT_Fx = data_MIT.Fx;

OUR_Fx = data_OUR.Fx;

%% Distance from true solution

MIT_DistSol = data_MIT.DistSol;

OUR_DistSol = data_OUR.DistSol;

%% Distance between methods
%% (identical for both methods)

DistBetween = data_MIT.DistBetween;

%% ============================================================
%% PLOT 1 : NUMBER OF FAILED SIMULATIONS
%% ============================================================

FailurePlot = figure();

failure_counts = [
    MIT_failures;
    OUR_failures
];

bar(failure_counts);

set(gca, ...
    'XTickLabel', {'MIT-LBFGSB', 'OUR-LBFGSB'});

ylabel('Number of failed simulations');

title('Non-converged simulations');

grid on;

saveas(FailurePlot, ...
       fullfile(folder, 'FailureHistogram.png'));

%% ============================================================
%% PLOT 2 : NUMBER OF ITERATIONS
%% ============================================================

NumIterPlot = figure();

NumIter_all = [
    MIT_NumIter;
    OUR_NumIter
];

group_iter = [
    repmat({'MIT-LBFGSB'}, length(MIT_NumIter), 1)
    repmat({'OUR-LBFGSB'}, length(OUR_NumIter), 1)
];

boxplot(NumIter_all, group_iter);

ylabel('Number of iterations');

title('Comparison of number of iterations');

grid on;

saveas(NumIterPlot, ...
       fullfile(folder, 'NumIterComparison.png'));

%% ============================================================
%% PLOT 3 : COMPUTATION TIME
%% ============================================================

CompTimePlot = figure();

CompTime_all = [
    MIT_CompTime;
    OUR_CompTime
];

group_time = [
    repmat({'MIT-LBFGSB'}, length(MIT_CompTime), 1)
    repmat({'OUR-LBFGSB'}, length(OUR_CompTime), 1)
];

boxplot(CompTime_all, group_time);

ylabel('Computation time [ms]');

title('Comparison of computation time');

grid on;

saveas(CompTimePlot, ...
       fullfile(folder, 'CompTimeComparison.png'));

%% ============================================================
%% PLOT 4 : DISTANCE FROM SOLUTION
%% ============================================================

DistSolPlot = figure();

DistSol_all = [
    MIT_DistSol;
    OUR_DistSol
];

group_distsol = [
    repmat({'MIT-LBFGSB'}, length(MIT_DistSol), 1)
    repmat({'OUR-LBFGSB'}, length(OUR_DistSol), 1)
];

boxplot(DistSol_all, group_distsol);

ylabel('Distance from correct solution');

title('Comparison of distance from correct solution');

grid on;

saveas(DistSolPlot, ...
       fullfile(folder, 'DistSolComparison.png'));

%% ============================================================
%% PLOT 5 : DISTANCE BETWEEN METHODS
%% ONLY KEEP TESTS WHERE BOTH METHODS CONVERGED
%% ============================================================

DistBetweenPlot = figure();

%% Extract original tables again
%% (before convergence filtering)

data_MIT_all = data(strcmp(data.OptId, 'MIT'), :);

data_OUR_all = data(strcmp(data.OptId, 'OUR'), :);

%% Logical masks for convergence

MIT_conv = data_MIT_all.DistSol <= tol;

OUR_conv = data_OUR_all.DistSol <= tol;

%% Keep only tests where BOTH methods converged

valid_tests = MIT_conv & OUR_conv;

%% Extract valid distances

DistBetween_valid = ...
    data_OUR_all.DistBetween(valid_tests);

group_between = ...
    repmat({'OUR-LBFGSB'}, length(DistBetween_valid), 1);

boxplot(DistBetween_valid, group_between);

ylabel('Distance from MIT solution');

title('Comparison of distance from MIT solution');

grid on;

saveas(DistBetweenPlot, ...
       fullfile(folder, 'DistBetweenComparison.png'));