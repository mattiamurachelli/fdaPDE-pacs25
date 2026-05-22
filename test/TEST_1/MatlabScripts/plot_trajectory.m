%% SCRIPT TO PLOT THE TRAJECTORIES OF THE MINIMIZATION
%% OF THE ROSENBROCK FUNCTION
%%
%% METHODS:
%% 1) OUR-LBFGSB
%% 2) OUR-ALM
%% 3) OUR-SQP

% CREATE PLOTS FOLDER
folder = '../Plots';

if ~exist(folder, 'dir')
    mkdir(folder)
end

%% READ TRAJECTORY FILES

traj_lbfgsb = readtable( ...
    '../trajectory_lbfgsb_rosenbrock.csv', ...
    'Delimiter', ',', ...
    'ReadVariableNames', true);

traj_alm = readtable( ...
    '../trajectory_alm_rosenbrock.csv', ...
    'Delimiter', ',', ...
    'ReadVariableNames', true);

traj_sqp = readtable( ...
    '../trajectory_sqp_rosenbrock.csv', ...
    'Delimiter', ',', ...
    'ReadVariableNames', true);

%% PARSE 2D POINTS (FORMAT: "x;y")

points_lbfgsb = parsePoints(traj_lbfgsb.Point);
points_alm    = parsePoints(traj_alm.Point);
points_sqp    = parsePoints(traj_sqp.Point);

%% EXTRACT X AND Y COMPONENTS

% OUR-LBFGSB
X_lbfgsb = points_lbfgsb(:,1);
Y_lbfgsb = points_lbfgsb(:,2);

% OUR-ALM
X_alm = points_alm(:,1);
Y_alm = points_alm(:,2);

% OUR-SQP
X_sqp = points_sqp(:,1);
Y_sqp = points_sqp(:,2);

%% DEFINE ROSENBROCK FUNCTION

rosenbrock = @(x,y) (1 - x).^2 + 100 * (y - x.^2).^2;

%% COMPUTE GLOBAL PLOT RANGE

all_x = [
    X_lbfgsb;
    X_alm;
    X_sqp
];

all_y = [
    Y_lbfgsb;
    Y_alm;
    Y_sqp
];

x_range = linspace(min(all_x)-1, max(all_x)+1, 200);
y_range = linspace(min(all_y)-1, max(all_y)+1, 200);

[X, Y] = meshgrid(x_range, y_range);

Z = rosenbrock(X, Y);

%% PLOT CONTOUR LINES

contourPlot = figure;

contour(X, Y, Z, logspace(-1,5,20), ...
        'LineWidth', 1.0);

hold on;

colormap('turbo');

colorbar;

xlabel('x');

ylabel('y');

title('Optimizer trajectories on Rosenbrock function');

grid on;

%% OVERLAY TRAJECTORIES

% OUR-LBFGSB
plot(X_lbfgsb, Y_lbfgsb, ...
     'r-o', ...
     'LineWidth', 1.5, ...
     'MarkerSize', 4, ...
     'MarkerFaceColor', 'r');

% OUR-ALM
plot(X_alm, Y_alm, ...
     'b-o', ...
     'LineWidth', 1.5, ...
     'MarkerSize', 4, ...
     'MarkerFaceColor', 'b');

% OUR-SQP
plot(X_sqp, Y_sqp, ...
     'g-o', ...
     'LineWidth', 1.5, ...
     'MarkerSize', 4, ...
     'MarkerFaceColor', 'g');

%% LEGEND

legend( ...
    'Contours', ...
    'OUR-LBFGSB', ...
    'OUR-ALM', ...
    'OUR-SQP', ...
    'Location', 'best');

%% SAVE FIGURE

saveas(contourPlot, ...
       fullfile(folder, "TrajectoryComparison.png"));