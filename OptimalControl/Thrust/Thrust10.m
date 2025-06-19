clc; clear; close all

%% Physical parameters

g = 9.81; % gravity [m/s^2]
rho = 1.225; % air density [kg/m^3]
Sw = 0.6; % aerodynamic surface [m^2]
mstruct = 4.5; % structural mass [kg]
Cd0 = 0.03; 
k = 0.05;
Clalpha = 0.1*360/(2*pi);  
Cl0 = 0;

%% Discretization

N = 200; % number of nodes
% load("resultsN300.mat")
% u_opt300 = u_opt; clear u_opt
% alpha0 = interp1(linspace(0,u_opt300(1),300),u_opt300(2:301),linspace(0,u_opt300(1),N));
% T0 = interp1(linspace(0,u_opt300(1),300),u_opt300(302:end),linspace(0,u_opt300(1),N));

%% Initial state

x1_0 = 0; x2_0 = 0; v0 = 15; gamma0 = deg2rad(30); m0 = 5;
t0 = 0; % Initial time [s]
tf0 = 10; % Initial guess for final time [s]
alpha0 = deg2rad(5); % Initial guess for the angle of attack [rad]
T0 = 8;
mfuel = m0 - mstruct; % fuel mass [kg]

%% Control vector

u0 = [tf0 alpha0*ones(1,N) T0*ones(1,N)];

%% Bounds for the control

alpha_min = deg2rad(-10);
alpha_max = deg2rad(10);
tf_min = 0.1;
tf_max = 15;
Tmin = 0;
Tmax = 10;

lb = [tf_min alpha_min*ones(1,N) Tmin*ones(1,N)]; % lower bound
ub = [tf_max alpha_max*ones(1,N) Tmax*ones(1,N)]; % upper bound

%% Cost function: maximize range

cost = @(u) f_cost(u, v0, gamma0, x1_0, x2_0, g, Cl0, Clalpha, Cd0, k, Sw, N, m0, rho);

%% Restrictions

nonlcon_fun = @(u) nonlcon(u, v0, x1_0, x2_0, gamma0, g, Cl0, Clalpha, Cd0, k, Sw, N, m0, rho, mstruct);

%% Optimization solver

options = optimoptions('fmincon','Algorithm', 'interior-point', ...
    'TolCon', 1e-9,'Display', 'iter','MaxFunctionEvaluations',1e6, ...
    'ScaleProblem',true,'StepTolerance',1e-6,'MaxIter',2000);

[u_opt, ~] = fmincon(cost, u0, [], [], [], [], lb, ub, nonlcon_fun, options);

%% Final results

tf = u_opt(1);
alpha = u_opt(2:N+1);
T = u_opt(N+2:end);
t_span = linspace(t0, tf, N);
dt = tf / (N-1);

z = trapezoidal(tf, alpha, N, x1_0, x2_0, v0, gamma0, g, Cl0, Clalpha, Cd0, k, Sw, m0, rho, T); % State vector solution
x = z(1,:); % Horizontal displacement [m]
y = z(2,:); % Vertical displacement [m]
v = z(3,:); % Velocity [m/s]
gamma = z(4,:); % Flight path angle [°]
m = z(5,:);

q = 0.5*rho.*v.^2;
Cl = Cl0+Clalpha.*alpha;
Cd = Cd0+k.*Cl.^2;
L = q.*Sw.*Cl;
D = q.*Sw.*Cd;

disp(["Maximum distance [m] = ", num2str(x(end))])
disp(["Final Time [s]: ", num2str(tf)])

[c_val, ceq_val] = nonlcon(u_opt, v0, x1_0, x2_0, gamma0, g, Cl0, Clalpha, Cd0, k, Sw, N, m0, rho, mstruct);
max_ineq = max(c_val);
final_eq = ceq_val;

%disp(['Max. violation of inequality constraint: ', num2str(max_ineq)])
%disp(['Violation of equality constraint y(tf)=0: ', num2str(final_eq)])

%% Postprocess

% Optimal trajectory 
figure
%plot(x, y, 'LineWidth', 5,'DisplayName','Trajectory','HandleVisibility','off') 
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
patch([x  NaN], [y  NaN], [T  NaN],'EdgeColor','interp', 'LineWidth', 4, ...
      'FaceColor','none','HandleVisibility','off');
colormap("turbo")
c = colorbar;
c.Label.String = 'Thrust [N]';

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
plot(t_span(1:end-1), rad2deg(gamma(1:end-1)), 'LineWidth', 2)
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

% Thrust vs. time

figure
plot(t_span,T,'DisplayName','Thrust','LineWidth',2)
grid on
xlabel('Time [s]','Interpreter','latex')
ylabel('Thrust [N]','Interpreter','latex')
title('Thrust profile','Interpreter','latex')

% Mass vs. time

figure
plot(t_span, m,'LineWidth',2)
grid on
xlabel('Time [s]','Interpreter','latex')
ylabel('Thrust [N]','Interpreter','latex')
title('Mass over time','Interpreter','latex')

%% ---------- FUNCTIONS ----------

function [J] = f_cost(u, v0, gamma0, x1_0, x2_0, g, Cl0, Clalpha, Cd0, k, Sw, N, m0, rho)
    tf = u(1);
    alpha = u(2:N+1);
    T = u(N+2:end);

    y = trapezoidal(tf, alpha, N, x1_0, x2_0, v0, gamma0, g, Cl0, Clalpha, Cd0, k, Sw, m0, rho, T);
    J = -y(1,end);  % maximize range
end

function [c, ceq] = nonlcon(u, v0, x1_0, x2_0, gamma0, g, Cl0, Clalpha, Cd0, k, Sw, N, m0, rho, mstruct)
    tf = u(1);
    alpha = u(2:N+1);
    T = u(N+2:end);
    y = trapezoidal(tf, alpha, N, x1_0, x2_0, v0, gamma0, g, Cl0, Clalpha, Cd0, k, Sw, m0, rho, T);
    c = [-y(2,:); mstruct-y(5,:)];         % y >= 0 
    ceq = [y(2,end); y(5,end)-mstruct];      % y(tf) = 0
end

function y = trapezoidal(tf, alpha, N, x1_0, x2_0, v0, gamma0, g, Cl0, Clalpha, Cd0, k, Sw, m0, rho, T)
    t_span = linspace(0, tf, N);
    dt = tf / (N-1);
    y = zeros(5,N);
    y(:,1) = [x1_0; x2_0; v0; gamma0; m0];
    for i = 1:N-1
        ti = t_span(i);
        alphai = alpha(i);
        alphai1 = alpha(i+1);
        Ti = T(i);
        Ti1 = T(i+1);
        yi = y(:,i);

        fi = dynamics_vector(ti, yi, alphai, g, Cl0, Clalpha, Cd0, k, Sw, rho, Ti);
        y_pred = yi + dt*fi;
        fi1 = dynamics_vector(t_span(i+1), y_pred, alphai1, g, Cl0, Clalpha, Cd0, k, Sw, rho, Ti1);
        y(:,i+1) = yi + (dt/2)*(fi + fi1);
    end
end

function dydt = dynamics_vector(~, y, alpha, g, Cl0, Clalpha, Cd0, k, Sw, rho, T)
    v = y(3);
    gamma = y(4);
    m = y(5);

    Cl =  Cl0+Clalpha*alpha;
    Cd =  Cd0+k*Cl^2;
    q = 0.5*rho*v^2;
    D =  q*Sw*Cd;
    L =  q*Sw*Cl;

    dydt = [
        v*cos(gamma);
        v*sin(gamma);
        -g*sin(gamma) - D/m + T/m;
        -(g/v)*cos(gamma) + L/(m*v);
        -T/20;
    ];
end
