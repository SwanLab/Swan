function main_constraint_drag
    persistent u_history J_history;
    clc; close all

    g = 9.81; % Gravity [m/s^2]
    aerodata.rho = 1.225; % Air density [kg/m^3]
    aerodata.Sw = 2; % Wing/Aerodynamic surface [m^2]
    aerodata.CD0 = 0.01; % Drag coefficient
    aerodata.m = 1; % Object mass [kg]

    N = 100; % Discretization

    x1_0 = 0; x2_0 = 0; v0 = 15;
    gamma0 = deg2rad(40);
    t0 = 0; tf = 30;
    alpha0 = deg2rad(10);
    u0 = [tf; alpha0];
    lb = [0 0.01]; % Lower bounds for the control
    ub = [50 deg2rad(10)]; % Upper bounds for the control

    u0 = [u0(1) ones(1,N)*u0(2)];
    lb = [lb(1) ones(1,N)*lb(2)];
    ub = [ub(1) ones(1,N)*ub(2)];

    cost = @(u) f_cost(u, g, t0, v0, gamma0, x1_0, x2_0, aerodata, N);
    constraint = @(u) groundConstraint(u, g, t0, v0, x1_0, x2_0, gamma0, aerodata, N);

    checkGradients(cost, u0);
    checkGradients(constraint, u0);

    options = optimoptions("fmincon", ...
        "OutputFcn", @store_fmincon, ...
        "Algorithm", "sqp", ...
        "DerivativeCheck", "on");

    [u_opt, ~] = fmincon(cost, u0, [], [], [], [], lb, ub, @(u) nonlcon(u,constraint), options);

    y0 = [x1_0 x2_0 v0 gamma0];
    t_span = linspace(t0, u_opt(1), 100);
    [~, y] = ode45(@(t, y) dynamics(t, y, tf, u_opt(2:N+1), g, aerodata, N), t_span, y0);

    disp(["Maximum distance [m] = ", num2str(y(end, 1))])
    disp(["Initial angle [°]: ", num2str(rad2deg(u_opt(2)))])
    disp(["Final Time [s]: ", num2str(u_opt(1))])

    figure; plot(y(:,1), y(:,2), 'b-', 'LineWidth', 2);
    xlabel("Horizontal distance [m]"); ylabel("Vertical distance [m]"); grid on;

    figure; plot(u_history(:,1), 'b-o', 'LineWidth', 2);
    xlabel('Iteration'); ylabel('Final time [s]'); grid on;

    figure; plot(rad2deg(u_history(:,2)), 'r-o', 'LineWidth', 2);
    xlabel('Iteration'); ylabel('Initial angle [°]'); grid on;

    figure; plot(J_history, 'o-', 'LineWidth', 2);
    xlabel('Iteration'); ylabel('Cost function J'); grid on;

    vars = whos;
    nodrag_results = struct();
    for i = 1:length(vars)
        nodrag_results.(vars(i).name) = eval(vars(i).name); 
    end
    save('nodrag_data.mat', 'nodrag_results');

    function stop = store_fmincon(u, optimValues, state)
        if optimValues.iteration == 0
            u_history = [];
            J_history = [];
        end
        if strcmp(state, 'iter')
            u_history = [u_history; u(:)'];
            J_history = [J_history; optimValues.fval];
            assignin('base', 'u_iterations', u_history);
            assignin('base', 'J_iterations', J_history);
        end
        stop = false;
    end
end

function [c, ceq, Dc, Dceq] = nonlcon(u,constraint)
    c = [];
    Dc = [];
    [ceq, Dceq] = constraint(u);
end

function [J, gradJ] = f_cost(u, g, t0, v0, gamma0, x1_0, x2_0, aerodata, N)
    y0 = [x1_0 x2_0 v0 gamma0];
    tf = u(1);
    alpha = u(2:N+1);
    t_span = linspace(t0, u(1), N);
    [t, y] = ode45(@(t, y) dynamics(t, y, tf, alpha, g, aerodata, N), t_span, y0);

    J = -y(end,1);
    dydt_final = dynamics(t_span(end), y(end,:)', tf, alpha, g, aerodata, N);

    pT = [1 0 0 0];
    [~, p] = ode45(@(t, p) p_ode(t, p, y, g, t_span, alpha, aerodata), flip(t_span), pT);
    q = p(end,:);

    v = interp1(t_span, y(:,3), t);
    DFdu = -0.1*aerodata.rho*v.^2*aerodata.Sw*pi^2/aerodata.m.*alpha';

    gradJ = [-dydt_final(1); DFdu*q(3)];
end

function [ceq, Dceq] = groundConstraint(u, g, t0, v0, x1_0, x2_0, gamma0, aerodata, N)
    y0 = [x1_0 x2_0 v0 gamma0]';
    tf = u(1);
    t_span = linspace(t0, tf, 100);
    alpha = u(2:N+1);
    [t, y] = ode45(@(t, y) dynamics(t, y, tf, alpha, g, aerodata, N), t_span, y0);

    ceq = y(end,2);
    dydt_final = dynamics(t_span(end), y(end,:)', tf, alpha, g, aerodata, N);

    pT = [0 -1 0 0];
    [~, p] = ode45(@(t, p) p_ode(t, p, y, g, t_span, alpha, aerodata), flip(t_span), pT);
    q = p(end,:);

    v = interp1(t_span, y(:,3), t);
    DFdu = -0.1*aerodata.rho*v.^2*aerodata.Sw*pi.^2./aerodata.m.*alpha';

    Dceq = [dydt_final(2); DFdu*q(3)];
end

function dpdt = p_ode(t, p, y, g, t_span, alpha, aerodata)
    v = interp1(t_span, y(:,3), t);
    gamma = interp1(t_span, y(:,4), t);
    alpha = interp1(t_span, alpha, t);
    
    CL = 2*pi*alpha;
    CD = aerodata.CD0 + 0.05*CL^2;

    J = [0 0 0 0;
         0 0 0 0;
         cos(gamma) sin(gamma) -aerodata.rho*aerodata.Sw*CD/aerodata.m*v (g/v^2)*cos(gamma);
         -v*sin(gamma) v*cos(gamma) -g*cos(gamma) (g/v)*sin(gamma)];
    dpdt = -J*p;
end

function dydt = dynamics(t, y, tf, alpha, g, aerodata, N)
    v = y(3);
    gamma = y(4);

    t_span = linspace(0, tf, N);
    alpha = interp1(t_span, alpha, t);
    CL = 2*pi*alpha;
    CD = aerodata.CD0 + 0.05*CL^2;

    dydt = [
        v*cos(gamma);
        v*sin(gamma);
        -g*sin(gamma)-0.5*aerodata.rho*CD*aerodata.Sw*v^2/aerodata.m;
        -(g/v)*cos(gamma)
    ];
end
