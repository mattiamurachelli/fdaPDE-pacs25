%% SCRIPT TO PLOT THE RESULTS OF THE SIMULATIONS OF
%% ALL OPTIMIZERS AGAINST MIT'S LBFGSB
%
% METHODS:
% 1) MIT-LBFGSB
% 2) OUR-LBFGSB
% 3) OUR-ALM
% 4) OUR-SQP
%

clear;
close all;
clc;

%% CREATE PLOTS FOLDER

folder = '../Plots';

if ~exist(folder, 'dir')
    mkdir(folder)
end

%% IMPORT CSV FILE

data = readtable('../simulations_rosenbrock.csv', ...
                 'ReadVariableNames', true, ...
                 'Delimiter', ',');

%% CONVERGENCE THRESHOLD

tol = 1e-3;

%% SEPARATE METHODS

data_MIT     = data(strcmp(data.OptId, 'MIT-LBFGSB'), :);
data_LBFGSB  = data(strcmp(data.OptId, 'OUR-LBFGSB'), :);
data_ALM     = data(strcmp(data.OptId, 'OUR-ALM'), :);
data_SQP     = data(strcmp(data.OptId, 'OUR-SQP'), :);

%% COUNT NON-CONVERGED SIMULATIONS

MIT_failures     = sum(data_MIT.DistSol     > tol);
LBFGSB_failures  = sum(data_LBFGSB.DistSol  > tol);
ALM_failures     = sum(data_ALM.DistSol     > tol);
SQP_failures     = sum(data_SQP.DistSol     > tol);

%% REMOVE NON-CONVERGED SIMULATIONS

data_MIT     = data_MIT(data_MIT.DistSol <= tol, :);
data_LBFGSB  = data_LBFGSB(data_LBFGSB.DistSol <= tol, :);
data_ALM     = data_ALM(data_ALM.DistSol <= tol, :);
data_SQP     = data_SQP(data_SQP.DistSol <= tol, :);

%% IMPORT DATA

% Number of iterations
MIT_NumIter     = data_MIT.NumIter;
LBFGSB_NumIter  = data_LBFGSB.NumIter;
ALM_NumIter     = data_ALM.NumIter;
SQP_NumIter     = data_SQP.NumIter;

% Computation time (convert ns -> ms)
MIT_CompTime     = data_MIT.CompTime    / 1e6;
LBFGSB_CompTime  = data_LBFGSB.CompTime / 1e6;
ALM_CompTime     = data_ALM.CompTime    / 1e6;
SQP_CompTime     = data_SQP.CompTime    / 1e6;

% Distance from solution
MIT_DistSol     = data_MIT.DistSol;
LBFGSB_DistSol  = data_LBFGSB.DistSol;
ALM_DistSol     = data_ALM.DistSol;
SQP_DistSol     = data_SQP.DistSol;

%% ============================================================
%% PLOT 1 : NON-CONVERGED SIMULATIONS
%% ============================================================

FailurePlot = figure();

failure_counts = [
    MIT_failures;
    LBFGSB_failures;
    ALM_failures;
    SQP_failures
];

bar(failure_counts);

set(gca, 'XTickLabel', ...
    {'MIT-LBFGSB', 'OUR-LBFGSB', 'OUR-ALM', 'OUR-SQP'});

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
    LBFGSB_NumIter;
    ALM_NumIter;
    SQP_NumIter
];

group_iter = [
    repmat({'MIT-LBFGSB'}, length(MIT_NumIter), 1)
    repmat({'OUR-LBFGSB'}, length(LBFGSB_NumIter), 1)
    repmat({'OUR-ALM'}, length(ALM_NumIter), 1)
    repmat({'OUR-SQP'}, length(SQP_NumIter), 1)
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
    LBFGSB_CompTime;
    ALM_CompTime;
    SQP_CompTime
];

group_time = [
    repmat({'MIT-LBFGSB'}, length(MIT_CompTime), 1)
    repmat({'OUR-LBFGSB'}, length(LBFGSB_CompTime), 1)
    repmat({'OUR-ALM'}, length(ALM_CompTime), 1)
    repmat({'OUR-SQP'}, length(SQP_CompTime), 1)
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
    LBFGSB_DistSol;
    ALM_DistSol;
    SQP_DistSol
];

group_distsol = [
    repmat({'MIT-LBFGSB'}, length(MIT_DistSol), 1)
    repmat({'OUR-LBFGSB'}, length(LBFGSB_DistSol), 1)
    repmat({'OUR-ALM'}, length(ALM_DistSol), 1)
    repmat({'OUR-SQP'}, length(SQP_DistSol), 1)
];

boxplot(DistSol_all, group_distsol);

ylabel('Distance from correct solution');

title('Comparison of distance from correct solution');

grid on;

saveas(DistSolPlot, ...
       fullfile(folder, 'DistSolComparison.png'));

%% ============================================================
%% PLOT 5 : DISTANCE BETWEEN METHODS AND MIT
%% ONLY KEEP TESTS WHERE BOTH METHODS CONVERGED
%% ============================================================

DistBetweenPlot = figure();

% Logical masks for convergence
MIT_conv     = data(strcmp(data.OptId, 'MIT-LBFGSB'), :).DistSol <= tol;
LBFGSB_conv  = data(strcmp(data.OptId, 'OUR-LBFGSB'), :).DistSol <= tol;
ALM_conv     = data(strcmp(data.OptId, 'OUR-ALM'), :).DistSol <= tol;
SQP_conv     = data(strcmp(data.OptId, 'OUR-SQP'), :).DistSol <= tol;

% Keep only tests where BOTH methods converged
valid_lbfgsb = MIT_conv & LBFGSB_conv;
valid_alm    = MIT_conv & ALM_conv;
valid_sqp    = MIT_conv & SQP_conv;

% Extract original tables again
data_MIT_all     = data(strcmp(data.OptId, 'MIT-LBFGSB'), :);
data_LBFGSB_all  = data(strcmp(data.OptId, 'OUR-LBFGSB'), :);
data_ALM_all     = data(strcmp(data.OptId, 'OUR-ALM'), :);
data_SQP_all     = data(strcmp(data.OptId, 'OUR-SQP'), :);

% Distances
DistBetween_LBFGSB = ...
    data_LBFGSB_all.DistBetween(valid_lbfgsb);

DistBetween_ALM = ...
    data_ALM_all.DistBetween(valid_alm);

DistBetween_SQP = ...
    data_SQP_all.DistBetween(valid_sqp);

DistBetween_all = [
    DistBetween_LBFGSB;
    DistBetween_ALM;
    DistBetween_SQP
];

group_between = [
    repmat({'OUR-LBFGSB'}, length(DistBetween_LBFGSB), 1)
    repmat({'OUR-ALM'}, length(DistBetween_ALM), 1)
    repmat({'OUR-SQP'}, length(DistBetween_SQP), 1)
];

boxplot(DistBetween_all, group_between);

ylabel('Distance from MIT solution');

title('Comparison of distance from MIT solution');

grid on;

saveas(DistBetweenPlot, ...
       fullfile(folder, 'DistBetweenComparison.png'));

