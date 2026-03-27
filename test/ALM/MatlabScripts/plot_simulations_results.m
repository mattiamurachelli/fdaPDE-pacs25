%% SCRIPT TO PLOT THE RESULTS OF THE SIMULATIONS OF OUR
%% ALM OPTIMIZER AGAINST MIT'S LBFGSB
% CREATE PLOTS FOLDER
folder = '../Plots';
if ~exist(folder, 'dir')
    mkdir(folder)
end
% IMPORT .CSV FILE
data = readtable('../simulations_rosenbrock.csv', 'ReadVariableNames', true, 'Delimiter', ',');
% SEPARATE by OptId( MIT vs OUR )
data_MIT = data(strcmp(data.OptId, 'MIT'), :);
data_OUR = data(strcmp(data.OptId, 'OUR'), :);
% IMPORT DATA
% NumIter
MIT_NumIter = data_MIT.NumIter;
OUR_NumIter = data_OUR.NumIter;
% CompTime (converting from ns to ms)
MIT_CompTime = data_MIT.CompTime/(1e6);
OUR_CompTime = data_OUR.CompTime/(1e6);
% MinPoint
MIT_MinPoint = parsePoints(data_MIT.MinPoint);
OUR_MinPoint = parsePoints(data_OUR.MinPoint);
% Fx
MIT_Fx = data_MIT.Fx;
OUR_Fx = data_OUR.Fx;
% DistSol
MIT_DistSol = data_MIT.DistSol;
OUR_DistSol = data_OUR.DistSol;
% DistBetween (it's identical for both, we only need one)
DistBetween = data_MIT.DistBetween;

% Plot NumIter
NumIterPlot = figure();
scatter(MIT_NumIter, OUR_NumIter, 'k*');
xlabel('MIT iterations');
ylabel('OUR iterations');
title('Comparison of iterations: MIT vs OUR');
grid on; 
hold on;
min_val = min([MIT_NumIter; OUR_NumIter]);
max_val = max([MIT_NumIter; OUR_NumIter]);
plot([min_val, max_val], [min_val, max_val], 'r--', 'LineWidth', 2);
legend('Data points', 'y = x', 'Location', 'best');
saveas(NumIterPlot, fullfile(folder, 'NumIterComparison.png'))

CompTimePlot = figure();
threshold = 1e2;
% Identify outliers
outliers = (MIT_CompTime > threshold) | (OUR_CompTime > threshold);
% Plot normal points
scatter(MIT_CompTime(~outliers), OUR_CompTime(~outliers), 'k*');
hold on;
% Plot outliers at the boundary
scatter(min(MIT_CompTime(outliers), threshold), ...
        min(OUR_CompTime(outliers), threshold), ...
        'ro', 'LineWidth', 1.5);
xlabel('MIT computation time');
ylabel('OUR computation time');
title('Comparison of computation time: MIT vs OUR');
grid on;
% Force axis limits
xlim([0 threshold]);
ylim([0 threshold]);
% Reference line
plot([0 threshold], [0 threshold], 'r--', 'LineWidth', 2);
legend('Data points', 'Outliers (>1e2)', 'y = x', 'Location', 'best');
saveas(CompTimePlot, fullfile(folder, 'CompTimeComparison.png'));

% Plot DistSolution
DistSolPlot = figure();
DistSol_all = [data_MIT.DistSol, data_OUR.DistSol];
boxplot(DistSol_all, 'Labels', {'MIT', 'OUR'});
ylabel('Distance from correct solution');
title('Comparison of distance from correct solution');
saveas(DistSolPlot, fullfile(folder, 'DistSolComparison.png'))

% Plot DistBetween
DistBetweenPlot = figure();
boxplot(DistBetween);
ylabel('Distance between the two solutions');
title('Comparison of distance between the two solutions');
saveas(DistBetweenPlot, fullfile(folder, 'DistBetweenComparison.png'));