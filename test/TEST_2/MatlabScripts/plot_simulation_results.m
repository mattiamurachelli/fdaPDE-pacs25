%% SCRIPT TO PLOT THE RESULTS OF THE SIMULATIONS OF ALM
%% AND SQP OPTIMIZERS ON PROBLEM HS071
% CREATE PLOTS FOLDER
folder = '../Plots';
if ~exist(folder, 'dir')
    mkdir(folder)
end
% IMPORT .CSV FILE
data = readtable('../simulations_hs071.csv', 'ReadVariableNames', true, 'Delimiter', ',');
% SEPARATE by OptId( ALM vs SQP )
data_ALM = data(strcmp(data.OptId, 'ALM'), :);
data_SQP = data(strcmp(data.OptId, 'SQP'), :);
% IMPORT DATA
% NumIter
ALM_NumIter = data_ALM.NumIter;
SQP_NumIter = data_SQP.NumIter;
% CompTime (converting from ns to ms)
ALM_CompTime = data_ALM.CompTime/(1e6);
SQP_CompTime = data_SQP.CompTime/(1e6);
% MinPoint
ALM_MinPoint = parsePoints(data_ALM.MinPoint);
SQP_MinPoint = parsePoints(data_SQP.MinPoint);
% Fx
ALM_Fx = data_ALM.Fx;
SQP_Fx = data_SQP.Fx;
% DistSol
ALM_DistSol = data_ALM.DistSol;
SQP_DistSol = data_SQP.DistSol;
% DistBetween (it's identical for both, we only need one)
DistBetween = data_ALM.DistBetween;

% Plot NumIter
NumIterPlot = figure();
NumIter_all = [
    ALM_NumIter;
    SQP_NumIter
];
group_iter = [
    repmat({'OUR-ALM'}, length(ALM_NumIter), 1)
    repmat({'OUR-SQP'}, length(SQP_NumIter), 1)
];
boxplot(NumIter_all, group_iter);
ylabel('Number of iterations');
title('Comparison of number of iterations');
grid on;

saveas(NumIterPlot, fullfile(folder, 'NumIterComparison.png'));

% Plot CompTime
CompTimePlot = figure();
CompTime_all = [
    ALM_CompTime;
    SQP_CompTime
];
group_time = [
    repmat({'OUR-ALM'}, length(ALM_CompTime), 1)
    repmat({'OUR-SQP'}, length(SQP_CompTime), 1)
];
boxplot(CompTime_all, group_time);
ylabel('Computation time [ms]');
title('Comparison of computation time');
grid on;
saveas(CompTimePlot, fullfile(folder, 'CompTimeComparison.png'));

% Plot DistSolution
DistSolPlot = figure();
DistSol_all = [data_ALM.DistSol, data_SQP.DistSol];
boxplot(DistSol_all, 'Labels', {'ALM', 'SQP'});
ylabel('Distance from correct solution');
title('Comparison of distance from correct solution');
saveas(DistSolPlot, fullfile(folder, 'DistSolComparison.png'))

% Plot DistBetween
DistBetweenPlot = figure();
boxplot(DistBetween);
ylabel('Distance between the two solutions');
title('Comparison of distance between the two solutions');
saveas(DistBetweenPlot, fullfile(folder, 'DistBetweenComparison.png'));