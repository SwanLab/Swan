% INDIRECT METHOD => Shooting technique to solve boundary problem 
%                    Solver for root finding -> fsolve


clc;
clear;

global error_x_history error_y_history
error_x_history = [];
error_y_history = [];

% Initial guess for [k, tf]

guess = [0, 40]; % k is c3 in the notes. tf is the final time

% Root-finding (using fsolve) to compute optimal values of k and tf 

options = optimoptions('fsolve', 'Display', 'iter', 'FunctionTolerance', 1e-6);

solution = fsolve(@residual_function, guess, options);

k_sol = solution(1);
tf_sol = solution(2);

fprintf('\nSolution found:\n');
fprintf('k = %.6f\n', k_sol);
fprintf('tf = %.6f\n', tf_sol);

% Integrate final trajectory once optimal k and tf are obtained

[x_sol, y_sol, t_sol] = integrate_trajectory(k_sol, tf_sol);

% Results plot

set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

figure;
plot(x_sol, y_sol, 'b-', 'LineWidth', 2);
hold on;
plot(0, 0, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g'); % Start
plot(1000, 200, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); % Target
grid on;
xlabel('x');
ylabel('y');
title('Optimal Trajectory Found Using Shooting Method');
legend('Trajectory', 'Start', 'Target');
ylim([0 250])
xlim([0 1000])

t_plot = linspace(0, tf_sol, 1000);
alpha_plot = atan(-0.5 * t_plot + k_sol);

figure;
plot(t_plot, rad2deg(alpha_plot), 'r', 'LineWidth', 2);
grid on;
xlabel('Time [s]');
ylabel('$\alpha$ [deg]');
title('Control Angle $\alpha$ Over Time');

V = 30;
xdot = V * cos(alpha_plot);
ydot = V * sin(alpha_plot);

figure;
plot(t_plot, xdot, 'b', 'LineWidth', 2);
hold on;
plot(t_plot, ydot, 'm--', 'LineWidth', 2);
grid on;
xlabel('Time [s]');
ylabel('Velocity Components [m/s]');
legend('$v_x$', '$v_y$');
title('Velocity Components Over Time');

figure;
subplot(2,1,1);
plot(1:length(error_x_history), error_x_history, 'b-o', 'LineWidth', 2);
grid on;
xlabel('Iteration');
ylabel('Error in x [m]');
title('Error in $x(t_f)$ Over FSOLVE Iterations');

subplot(2,1,2);
plot(1:length(error_y_history), error_y_history, 'r-s', 'LineWidth', 2);
grid on;
xlabel('Iteration');
ylabel('Error in y [m]');
title('Error in $y(t_f)$ Over FSOLVE Iterations');


% Function Definitions

function F = residual_function(params)
    global error_x_history error_y_history
    k = params(1);
    tf = params(2);
    
    % Integrate system from t=0 to t=tf
    [x_sol, y_sol, ~] = integrate_trajectory(k, tf);
    
    x_tf = x_sol(end);
    y_tf = y_sol(end);
    
    % Errores con respecto al objetivo
    ex = x_tf - 1000;
    ey = y_tf - 200;

    % Guardar los errores por separado
    error_x_history(end+1) = ex;
    error_y_history(end+1) = ey;

    % Boundary errors
    F = [x_tf - 1000; y_tf - 200];
end

function [x_vals, y_vals, t_vals] = integrate_trajectory(k, tf)
    V = 30;
    z0 = [0; 0]; % Initial conditions: x(0)=0, y(0)=0
    
    % Define ODE system
    odefun = @(t, z) [
        V * cos(atan(-0.5 * t + k)) + 0.5 * z(2); % xdot
        V * sin(atan(-0.5 * t + k))               % ydot
    ];
    
    % Integrate with ode45
    [t_vals, z] = ode45(odefun, [0 tf], z0);
    x_vals = z(:,1);
    y_vals = z(:,2);
end
