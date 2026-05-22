%% SCRIPT TO PLOT THE RESULTS OF THE SIMULATIONS OF ALM
%% VS MIT LBFGSB OPTIMIZER
% CREATE PLOTS FOLDER
folder = '../Plots';
if ~exist(folder, 'dir')
    mkdir(folder)
end
% IMPORT .CSV FILE
data = readtable('../simulations_rosenbrock.csv', 'ReadVariableNames', true, 'Delimiter', ',');
% SEPARATE by OptId( MIT vs ALM )
data_MIT = data(strcmp(data.OptId, 'MIT'), :);
data_ALM = data(strcmp(data.OptId, 'ALM'), :);
% IMPORT DATA
% NumIter
MIT_NumIter = data_MIT.NumIter;
ALM_NumIter = data_ALM.NumIter;
% CompTime (converting from ns to ms)
MIT_CompTime = data_MIT.CompTime/(1e6);
ALM_CompTime = data_ALM.CompTime/(1e6);
% MinPoint
MIT_MinPoint = parsePoints(data_MIT.MinPoint);
ALM_MinPoint = parsePoints(data_ALM.MinPoint);
% Fx
MIT_Fx = data_MIT.Fx;
ALM_Fx = data_ALM.Fx;
% DistSol
MIT_DistSol = data_MIT.DistSol;
ALM_DistSol = data_ALM.DistSol;
% DistBetween (it's identical for both, we only need one)
DistBetween = data_MIT.DistBetween;

% Plot NumIter
NumIterPlot = figure();
NumIter_all = [
    MIT_NumIter;
    ALM_NumIter
];
group_iter = [
    repmat({'OUR-MIT'}, length(MIT_NumIter), 1)
    repmat({'OUR-ALM'}, length(ALM_NumIter), 1)
];
boxplot(NumIter_all, group_iter);
ylabel('Number of iterations');
title('Comparison of number of iterations');
grid on;

saveas(NumIterPlot, fullfile(folder, 'NumIterComparison.png'));

% Plot CompTime
CompTimePlot = figure();
CompTime_all = [
    MIT_CompTime;
    ALM_CompTime
];
group_time = [
    repmat({'OUR-MIT'}, length(MIT_CompTime), 1)
    repmat({'OUR-ALM'}, length(ALM_CompTime), 1)
];
boxplot(CompTime_all, group_time);
ylabel('Computation time [ms]');
title('Comparison of computation time');
grid on;
saveas(CompTimePlot, fullfile(folder, 'CompTimeComparison.png'));

% Plot DistSolution
DistSolPlot = figure();
DistSol_all = [data_MIT.DistSol, data_ALM.DistSol];
boxplot(DistSol_all, 'Labels', {'MIT', 'ALM'});
ylabel('Distance from correct solution');
title('Comparison of distance from correct solution');
saveas(DistSolPlot, fullfile(folder, 'DistSolComparison.png'))

% Plot DistBetween
DistBetweenPlot = figure();
boxplot(DistBetween);
ylabel('Distance between the two solutions');
title('Comparison of distance between the two solutions');
saveas(DistBetweenPlot, fullfile(folder, 'DistBetweenComparison.png'));