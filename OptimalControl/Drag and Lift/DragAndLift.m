% Launch with drag and lift

clc; clear; close all

%% Physical parameters

g = 9.81;
m = 5.0;
rho = 1.225;
S = 0.6;
Cd0 = 0.03;
k = 0.05;
Cl0 = 0.0;
Clalpha = 0.1*360/(2*pi);

%% Discretization

N = 50;             % Number of nodes

%% Initial state

x0 = 0; y0 = 0; v0 = 15; gamma0 = deg2rad(30); % State initial conditions
alpha0 = deg2rad(5); % Initial guess for angle of attack [rad]  
tf0 = 7;            % Initial guess for final time [s]

%% Decision variables vector

z0 = repmat([x0; y0; v0; gamma0; alpha0], N, 1);
z0 = [z0; tf0];

%% Bounds for the variables

alpha_min = deg2rad(-10);
alpha_max = deg2rad(10);
tf_min = 2;
tf_max = 30;

lb = repmat([-Inf; 0; 0; -pi/2; alpha_min], N, 1); % Lower bounds
ub = repmat([ Inf; Inf; Inf;  pi/2; alpha_max], N, 1); % Upper bounds
lb = [lb; tf_min];
ub = [ub; tf_max];

%% Cost function: maximize range

obj = @(z) -z((N-1)*5 + 1);  % Cost function x(t_f)

%% Restrictions

nonlcon = @(z) collocation_constraints_tf(z, N, rho, S, Cd0, k, Cl0, Clalpha, m, g, v0, gamma0);

%% Optimization

opts = optimoptions('fmincon', ...
    'Display','iter-detailed', ...
    'MaxFunctionEvaluations', 1e6, ...
    'MaxIterations', 1000, ...
    'StepTolerance', 1e-6, ...
    'ScaleProblem', true, ...
    'Algorithm', 'interior-point');

[z_opt, fval] = fmincon(obj, z0, [], [], [], [], lb, ub, nonlcon, opts);

%% Final results

tf_opt = z_opt(end); % Optimal final time [s]
Z = reshape(z_opt(1:end-1), 5, N);
x = Z(1,:); y = Z(2,:); v = Z(3,:); gamma = Z(4,:); alpha = Z(5,:);
t = linspace(0, tf_opt, N);

q = 0.5*rho.*v.^2;
Cl = Cl0+Clalpha.*alpha;
Cd = Cd0+k.*Cl.^2;
L = q.*S.*Cl;
D = q.*S.*Cd;

fprintf('Máximum range: %.2f m\n', x(end));
fprintf('Optimal final time: %.2f s\n', tf_opt);

%% Postprocess

set(groot,'defaultLineLineWidth',1.5);   % thicker lines for all plots

% Optimal trajectory 
figure
plot(x, y, 'LineWidth', 2,'DisplayName','Trajectory') 
hold on
plot(x(1), y(1), 'go', 'MarkerSize', 7, 'MarkerFaceColor', 'g', 'DisplayName','Launch');
plot(x(end),y(end),'ro', 'MarkerSize', 7, 'MarkerFaceColor','r','DisplayName','Impact');
xlim([0 1.05*max(x)])
ylim([min(0,1.05*min(y))  1.05*max(y)])
grid on; box on
xlabel('x [m]','Interpreter','latex')
ylabel('y [m]','Interpreter','latex')
title('Optimal trajectory','Interpreter','latex')
legend('Location','best')
set(gca,'FontSize',12)

% Angle-of-attack profile 
figure
stairs(t, rad2deg(alpha), 'Marker','o')
grid on; box on
xlim([0 t(end)])
ylim([min(rad2deg(alpha))-2  max(rad2deg(alpha))+2])
xlabel('Time [s]','Interpreter','latex')
ylabel('$\alpha$ [$^\circ$]','Interpreter','latex')
title('Angle of attack profile','Interpreter','latex')
set(gca,'FontSize',12)

% Speed vs. time
figure
plot(t, v, 'LineWidth', 2)
grid on; 
xlabel('Time [s]','Interpreter','latex')
ylabel('$v$ [m/s]','Interpreter','latex')
title('Velocity profile','Interpreter','latex')
set(gca,'FontSize',12)

% Flight-path angle vs. time 
figure
plot(t, rad2deg(gamma), 'LineWidth', 2)
grid on; box on
xlabel('Time [s]','Interpreter','latex')
ylabel('$\gamma$ [$^\circ$]','Interpreter','latex')
title('Flight path angle','Interpreter','latex')
set(gca,'FontSize',12)

% Lift & Drag vs. time 

figure
plot(t, L,'DisplayName','Lift','LineWidth',2); 
hold on
plot(t, D, 'DisplayName','Drag','LineWidth',2);
grid on; box on
xlabel('Time [s]','Interpreter','latex')
ylabel('Force [N]','Interpreter','latex')
title('Aerodynamic forces','Interpreter','latex')
legend('Location','best')
set(gca,'FontSize',12)

%% ---------- FUNCTIONS ----------

function [c, ceq] = collocation_constraints_tf(z, N, rho, S, Cd0, k, Cl0, Clalpha, m, g, v0, gamma0)
    tf = z(end);
    dt = tf / (N - 1);
    Z = reshape(z(1:end-1), 5, N);
    x = Z(1,:); y = Z(2,:); v = Z(3,:); gamma = Z(4,:); alpha = Z(5,:);

    ceq = [];
    for i = 1:N-1
        % States at i and i+1
        xi = x(i);     yi = y(i);     vi = v(i);     gi = gamma(i);     ai = alpha(i);
        xj = x(i+1);   yj = y(i+1);   vj = v(i+1);   gj = gamma(i+1);   aj = alpha(i+1);

        % Derivatives
        [dxi, dyi, dvi, dgi] = dynamics(xi, yi, vi, gi, ai, rho, S, Cd0, k, Cl0, Clalpha, m, g);
        [dxj, dyj, dvj, dgj] = dynamics(xj, yj, vj, gj, aj, rho, S, Cd0, k, Cl0, Clalpha, m, g);

        % Collocation equations (trapezoidal)
        ceq = [ceq;
            xj-xi-dt/2*(dxi + dxj);
            yj-yi-dt/2*(dyi + dyj);
            vj-vi-dt/2*(dvi + dvj);
            gj-gi-dt/2*(dgi + dgj)];
    end

    % Initial conditions
    ceq = [ceq;
        x(1)-0;
        y(1)-0;
        v(1)-v0;
        gamma(1)-gamma0];

    % Equality constraint: y(tf) = 0
    ceq = [ceq;
        y(end)];

    % Inequality constraint: y_i ≥ 0
    c = -y;

    % Optimal trajectory 
figure(100)
plot(x, y, 'LineWidth', 2,'DisplayName','Trajectory') 
hold on
plot(x(1), y(1), 'go', 'MarkerSize', 7, 'MarkerFaceColor', 'g', 'DisplayName','Launch');
plot(x(end),y(end),'ro', 'MarkerSize', 7, 'MarkerFaceColor','r','DisplayName','Impact');
xlim([0 1.05*max(x)+0.01])
ylim([min(0,1.05*min(y))  1.05*max(y)])
grid on; box on
xlabel('x [m]','Interpreter','latex')
ylabel('y [m]','Interpreter','latex')
title('Optimal trajectory','Interpreter','latex')
legend('Location','best')
set(gca,'FontSize',12)
hold off

end

function [dx, dy, dv, dgamma] = dynamics(x, y, v, gamma, alpha, rho, S, Cd0, k, Cl0, Clalpha, m, g)
    Cl = Cl0+Clalpha*alpha;
    Cd = Cd0+k* Cl^2;
    q = 0.5*rho*v^2;
    L = q*S*Cl;
    D = q*S*Cd;

    dx = v*cos(gamma);
    dy = v*sin(gamma);
    dv = -D/m-g*sin(gamma);
    dgamma = L/(m*v)-g*cos(gamma)/v;

end
