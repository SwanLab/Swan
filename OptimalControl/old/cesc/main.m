% DIRECT COLLOCATION METHOD => Transcription done by Falcon
%                              Trapezoidal discretization 
%                              NLP solver -> IPOPT

clear;
clc;

falcon.init();

%% Define States Controls and Parameter
x_vec = [...
    falcon.State('x',       -1000,    2000, 0.001);...
    falcon.State('y',       -1000,    2000, 0.001)];

u_vec = [...
    falcon.Control('u'  ,   -pi,    pi, 1)];
    
tf = falcon.Parameter('FinalTime', 20, 0, 1000, 0.01);

%% Define Optimal Control Problem
% Create new Problem Instance (Main Instance)
problem = falcon.Problem('Zermelo_wind');

% Specify Discretization
tau = linspace(0,1,201);

% Add a new Phase
phase = problem.addNewPhase(@source_zermelo, x_vec, tau, 0, tf);
phase.addNewControlGrid(u_vec, tau);


% Set Boundary Condition
phase.setInitialBoundaries([0,0]);
phase.setFinalBoundaries([1000,200]);

  
% Add Cost Function
problem.addNewLinearPointCost(tf);


% Prepare problem for solving
problem.Bake();

% Solve problem
solver = falcon.solver.ipopt(problem);
solver.Options.MajorIterLimit = 500;
solver.Options.MajorFeasTol   = 1e-6;
solver.Options.MajorOptTol    = 1e-6;

solver.Solve();

%% Plots

results = struct();
results.t = problem.RealTime;
results.tfinal = problem.Parameters.byName().FinalTime.Value;
results.x = x_vec.extractValuesStructFrom(problem);
results.u = u_vec.extractValuesStructFrom(problem);

figure();

ax = subplot(2, 2, 1);
plot(ax, results.t, rad2deg(results.u.u));
xlabel(ax, "t");
ylabel(ax, "alfa (deg)");

ax = subplot(2, 2, 2);
plot(ax, results.x.x, results.x.y);
xlabel(ax, "x");
ylabel(ax, "y");

% Plot with wind vector field

% Define the range for x and y
x = linspace(0, 1100, 20);
y = linspace(0, 250, 20);

% Create meshgrid for x and y
[X, Y] = meshgrid(x, y);

% Define the vector field
vx = 0.05 * Y;
vy= 0 * X;

set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

% Plot the vector field and the trajectory using quiver
figure;
quiver(X, Y, vx, vy, 'AutoScaleFactor', 0.75);
hold on; % Keep the figure for further plotting

% Plot the trajectory on top of the vector field
plot(results.x.x, results.x.y, 'b', 'LineWidth', 2);

% Plot starting point
plot(results.x.x(1), results.x.y(1), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');

% Plot target point
plot(results.x.x(end), results.x.y(end), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

% Add annotations and grid
xlim([0 1000]); % Setting x-axis limits
ylim([0 250]); % Setting y-axis limits
title('Wind Vector Field with Trajectory');
xlabel('X Coordinate');
ylabel('Y Coordinate');
legend('Wind Vectors', 'Trajectory', 'Location', 'best');
grid on;
hold off; % Release the figure

% Control plot

figure;
plot(results.t, rad2deg(results.u.u), 'r', 'LineWidth', 2);
grid on;
xlabel('Time [s]');
ylabel('\alpha [^\circ]');
title('Control Angle \alpha Over Time');

results.tfinal
