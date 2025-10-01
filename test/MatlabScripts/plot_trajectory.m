%% SCRIPT TO PLOT THE TRAJECTORY OF THE MINIMIZATION OF THE ROSENBROCK FUNCTION
%% OF OUR LBFGSB IMPLEMENTATION

% Read trajectory file
traj = readtable('trajectory_rosenbrock.csv', 'Delimiter', ',', 'ReadVariableNames', true);

% Parse 2D points (in "x;y" format)
points = parsePoints(traj.Point);

% Extract X and Y
X_traj = points(:,1);
Y_traj = points(:,2);

% Define Rosenbrock function
rosenbrock = @(x,y) (1 - x).^2 + 100*(y - x.^2).^2;

% Make a grid for contours
x_range = linspace(min(X_traj)-1, max(X_traj)+1, 200);
y_range = linspace(min(Y_traj)-1, max(Y_traj)+1, 200);
[X, Y] = meshgrid(x_range, y_range);
Z = rosenbrock(X,Y);

% Plot contour lines
contourPlot = figure;
contour(X, Y, Z, logspace(-1,5,20), 'LineWidth', 1.0); % logarithmic levels for better visualization
hold on;
colormap('turbo');
colorbar;
xlabel('x');
ylabel('y');
title('Optimizer trajectory on Rosenbrock function');
grid on;

% Overlay trajectory
plot(X_traj, Y_traj, 'r-o', 'LineWidth', 1, 'MarkerSize', 5, 'MarkerFaceColor', 'r');

% Highlight start and end
plot(X_traj(1), Y_traj(1), 'gs', 'MarkerSize', 8, 'MarkerFaceColor', 'g'); % start point
plot(X_traj(end), Y_traj(end), 'bs', 'MarkerSize', 8, 'MarkerFaceColor', 'b'); % final point

legend('Contours','Trajectory','Start','End', 'Location', 'best');
saveas(contourPlot, "Trajectory.png")
