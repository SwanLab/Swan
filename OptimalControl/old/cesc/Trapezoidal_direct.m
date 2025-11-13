% DIRECT COLLOCATION METHOD => Solved with own transcription
%                              Trapezoidal discretization
%                              NLP solver -> fmincon

clear;
clc;

% Problem inputs

N = 50;                         % number of mesh intervals
x0 = 0;  y0 = 0;                % initial state
xf = 1000;  yf = 200;           % target state
V = 30;                         % Velocity module

% Initial guesses

tf_guess = 40;                  % initial guess for final time
x_guess     = linspace(x0, xf, N+1);  % straight line
y_guess     = linspace(y0, yf, N+1);  % straight line
alpha_guess = atan2(diff(y_guess), diff(x_guess));
alpha_guess = [alpha_guess, alpha_guess(end)];   %  N+1

z0          = [x_guess, y_guess, alpha_guess, tf_guess]';  % column

% bounds for variables

lb = [-inf*ones(3*(N+1),1);  1];          % tf ≥ 1 s
ub = [ inf*ones(2*(N+1),1);  pi*ones(N+1,1);  1e4];  % tf ≤ 10 000 s

% fmincon solver options

opts = optimoptions('fmincon','Display','iter','MaxFunEvals',1e6,...
                    'MaxIter',1e3,'Algorithm','interior-point');

% solve NLP problem using fmincon

z_opt = fmincon(@objective, z0, [], [], [], [], ...
                lb, ub, @(z) constraints(z, xf, yf,V), opts);

% plot result

np1   = (length(z_opt)-1)/3;
x     = z_opt(1:np1);
y     = z_opt(np1+1:2*np1);
alpha = z_opt(2*np1+1:3*np1);
tf    = z_opt(end);
t     = linspace(0, tf, np1);

set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

figure;
plot(x, y, 'b.-');
xlabel('x'); ylabel('y'); grid on;
title('Optimal path');

alpha_deg = rad2deg(wrapToPi(alpha));

figure;
plot(t, alpha_deg, 'r.-','LineWidth', 2, 'MarkerSize', 10);
xlabel('Time [s]'); 
ylabel('$\alpha$ $[^\circ]$'); 
grid on;
title('Optimal heading angle');

% --- Trajectory plot ---
figure;
plot(x, y, 'b.-', 'LineWidth', 2, 'MarkerSize', 10); 
hold on;
plot(x(1), y(1), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g'); % Start
plot(x(end), y(end), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); % Target
grid on;
xlabel('x [m]');
ylabel('y [m]');
ylim([0 250])
xlim([0 1000])
title('Optimal Trajectory using Direct Collocation');
legend('Trajectory', 'Start', 'Target', 'Location', 'best');

% cost function (minimize final time)

function J = objective(z)
    J = z(end);
end

% constraints 
function [c, ceq] = constraints(z, xf, yf,V)

    % decode decision vector
    np1 = (length(z)-1)/3;     % np1 = N + 1
    N   = np1 - 1;             % number of intervals

    x     = z(1:np1);
    y     = z(np1+1:2*np1);
    alpha = z(2*np1+1:3*np1);
    tf    = z(end);
    h     = tf / N;

    % trapezoidal dynamics
    ceq_dyn = zeros(2*N,1);
    for i = 1:N
        f1_i   = V*cos(alpha(i))     + 0.5*y(i);
        f1_ip1 = V*cos(alpha(i+1))   + 0.5*y(i+1);

        f2_i   = V*sin(alpha(i));
        f2_ip1 = V*sin(alpha(i+1));

        ceq_dyn(i)     = x(i+1) - x(i) - h/2*(f1_i + f1_ip1);
        ceq_dyn(N+i)   = y(i+1) - y(i) - h/2*(f2_i + f2_ip1);
    end

    %  boundary conditions
    ceq_bc = [x(1); y(1); x(end)-xf; y(end)-yf];

    % Put together all the constraints 
    ceq = [ceq_dyn; ceq_bc];    % equality constraints
    c   = [];                   % no inequality constraints
end
