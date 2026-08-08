clc; clear; close all

%% Physical parameters

g = 9.81; % gravity [m/s^2]
rho = 1.225; % air density [kg/m^3]
Sw = 0.6; % aerodynamic surface [m^2]
Cd0 = 0.03; 
m = 5; % mass [kg]
k = 0.05;
Clalpha = 0.1*360/(2*pi);  
Cl0 = 0;

%% Discretization

N = 250; % number of nodes

%% Initial state

x1_0 = 0; x2_0 = 0; v0 = 15; gamma0 = deg2rad(30);
t0 = 0; % Initial time [s]
tf0 = 7; % Initial guess for final time [s]
alpha0 = deg2rad(5); % Initial guess for the angle of attack [rad]

%% Control vector

u0 = [tf0 alpha0*ones(1,N)]; 

%% Bounds for the control

alpha_min = deg2rad(-10);
alpha_max = deg2rad(10);
tf_min = 0.1;
tf_max = 15;

lb = [tf_min alpha_min*ones(1,N)]; % lower bound
ub = [tf_max alpha_max*ones(1,N)]; % upper bound

%% Cost function: maximize range

cost = @(u) f_cost(u, v0, gamma0, x1_0, x2_0, g, Cl0, Clalpha, Cd0, k, Sw, m, N, rho);

%% Restrictions

nonlcon_fun = @(u) nonlcon(u, v0, x1_0, x2_0, gamma0, g, Cl0, Clalpha, Cd0, k, Sw, m, N, rho);

%% Optimization solver

options = optimoptions('fmincon','Algorithm', 'interior-point', ...
    'TolCon', 1e-8,'Display', 'iter','MaxFunctionEvaluations',1e6, ...
    'ScaleProblem',true,'StepTolerance',1e-6,'MaxIter',2000);

[u_opt, ~] = fmincon(cost, u0, [], [], [], [], lb, ub, nonlcon_fun, options);

%% Final results

tf = u_opt(1);
alpha = u_opt(2:N+1);
t_span = linspace(t0, tf, N);
dt = tf / (N-1);

z = trapezoidal(tf, alpha, N, x1_0, x2_0, v0, gamma0, g, Cl0, Clalpha, Cd0, k, Sw, m, rho); % State vector solution
x = z(1,:); % Horizontal displacement [m]
y = z(2,:); % Vertical displacement [m]
v = z(3,:); % Velocity [m/s]
gamma = z(4,:); % Flight path angle [°]

q = 0.5*rho.*v.^2;
Cl = Cl0+Clalpha.*alpha;
Cd = Cd0+k.*Cl.^2;
L = q.*Sw.*Cl;
D = q.*Sw.*Cd;

disp(["Maximum distance [m] = ", num2str(x(end))])
disp(["Final Time [s]: ", num2str(tf)])

[c_val, ceq_val] = nonlcon(u_opt, v0, x1_0, x2_0, gamma0, g, Cl0, Clalpha, Cd0, k, Sw, m, N, rho);
max_ineq = max(c_val);
final_eq = ceq_val;

disp(['Max. violation of inequality constraint: ', num2str(max_ineq)])
disp(['Violation of equality constraint y(tf)=0: ', num2str(final_eq)])

validateGradient(u0, v0, gamma0, x1_0, x2_0, g, Cl0, Clalpha, Cd0, k, Sw, m, N, rho)

%% Postprocess

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
stairs(t_span, rad2deg(alpha), 'Marker','o')
grid on; box on
xlim([0 t_span(end-1)])
ylim([min(rad2deg(alpha))-2  max(rad2deg(alpha))+2])
xlabel('Time [s]','Interpreter','latex')
ylabel('$\alpha$ [$^\circ$]','Interpreter','latex')
title('Angle of attack profile','Interpreter','latex')
set(gca,'FontSize',12)

% Speed vs. time
figure
plot(t_span, v, 'LineWidth', 2)
grid on; 
xlabel('Time [s]','Interpreter','latex')
ylabel('$v$ [m/s]','Interpreter','latex')
title('Velocity profile','Interpreter','latex')
set(gca,'FontSize',12)

% Flight-path angle vs. time 
figure
plot(t_span, rad2deg(gamma), 'LineWidth', 2)
grid on; box on
xlabel('Time [s]','Interpreter','latex')
ylabel('$\gamma$ [$^\circ$]','Interpreter','latex')
title('Flight path angle','Interpreter','latex')
set(gca,'FontSize',12)

% Lift & Drag vs. time 

figure
plot(t_span(1:end-1), L(1:end-1),'DisplayName','Lift','LineWidth',2); 
hold on
plot(t_span(1:end-1), D(1:end-1), 'DisplayName','Drag','LineWidth',2);
grid on; box on
xlabel('Time [s]','Interpreter','latex')
ylabel('Force [N]','Interpreter','latex')
title('Aerodynamic forces','Interpreter','latex')
legend('Location','best')
set(gca,'FontSize',12)

%% ---------- FUNCTIONS ----------

function [J] = f_cost(u, v0, gamma0, x1_0, x2_0, g, Cl0, Clalpha, Cd0, k, Sw, m, N, rho)
    tf = u(1);
    alpha = u(2:N+1);

    y = trapezoidal(tf, alpha, N, x1_0, x2_0, v0, gamma0, g, Cl0, Clalpha, Cd0, k, Sw, m, rho);
    J = -y(1,end);  % maximize range
end

function [c, ceq] = nonlcon(u, v0, x1_0, x2_0, gamma0, g, Cl0, Clalpha, Cd0, k, Sw, m, N, rho)
    tf = u(1);
    alpha = u(2:N+1);
    y = trapezoidal(tf, alpha, N, x1_0, x2_0, v0, gamma0, g, Cl0, Clalpha, Cd0, k, Sw, m, rho);
    c = -y(2,:);         % y >= 0 
    ceq = y(2,end);      % y(tf) = 0
end

function y = trapezoidal(tf, alpha, N, x1_0, x2_0, v0, gamma0, g, Cl0, Clalpha, Cd0, k, Sw, m, rho)
    t_span = linspace(0, tf, N);
    dt = tf / (N-1);
    y = zeros(4,N);
    y(:,1) = [x1_0; x2_0; v0; gamma0];
    for i = 1:N-1
        ti = t_span(i);
        alphai = alpha(i);
        alphai1 = alpha(i+1);
        yi = y(:,i);

        fi = dynamics_vector(ti, yi, alphai, g, Cl0, Clalpha, Cd0, k, Sw, m, rho);
        y_pred = yi + dt*fi;
        fi1 = dynamics_vector(t_span(i+1), y_pred, alphai1, g, Cl0, Clalpha, Cd0, k, Sw, m, rho);
        y(:,i+1) = yi + (dt/2)*(fi + fi1);
    end
end

function dydt = dynamics_vector(~, y, alpha, g, Cl0, Clalpha, Cd0, k, Sw, m, rho)
    v = y(3);
    gamma = y(4);

    Cl =  Cl0+Clalpha*alpha;
    Cd =  Cd0+k*Cl^2;
    q = 0.5*rho*v^2;
    D =  q*Sw*Cd;
    L =  q*Sw*Cl;

    dydt = [
        v*cos(gamma);
        v*sin(gamma);
        -g*sin(gamma) - D/m;
        -(g/v)*cos(gamma) + L/(m*v)
    ];
end
