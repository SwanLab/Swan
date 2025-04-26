function main_constraint_nodrag
    persistent u_history J_history;
    clc; close all

    g = 9.81;
    lb = [0 0.01];
    ub = [10 pi/2];

    x1_0 = 0; x2_0 = 0; v0 = 15;
    gamma0 = deg2rad(40);
    t0 = 0; tf = 2;
    u0 = [tf; gamma0];

    cost = @(u) f_cost(u, g, t0, v0, x1_0, x2_0);
    constraint = @(u) groundConstraint(u, g, t0, v0, x1_0, x2_0);

    checkGradients(cost, u0);
    checkGradients(constraint, u0);

    options = optimoptions("fmincon", ...
        "OutputFcn", @store_fmincon, ...
        "Algorithm", "sqp", ...
        "DerivativeCheck", "on");

    [u_opt, ~] = fmincon(cost, u0, [], [], [], [], lb, ub, @(u) nonlcon(u,constraint), options);

    y0 = [x1_0 x2_0 v0 u_opt(2)];
    t_span = linspace(t0, u_opt(1), 100);
    [~, y] = ode45(@(t, y) dynamics(y, g), t_span, y0);

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

function [J, gradJ] = f_cost(u, g, t0, v0, x1_0, x2_0)
    y0 = [x1_0 x2_0 v0 u(2)];
    t_span = linspace(t0, u(1), 1000);
    [~, y] = ode45(@(t, y) dynamics(y, g), t_span, y0);

    J = -y(end,1);
    dydt_final = dynamics(y(end,:)', g);

    pT = [1 0 0 0];
    [~, p] = ode45(@(t, p) p_ode(t, p, y, g, t_span), flip(t_span), pT);
    q = p(end,:);

    gradJ = [-dydt_final(1); -q(4)];
end

function [ceq, Dceq] = groundConstraint(u, g, t0, v0, x1_0, x2_0)
    y0 = [x1_0 x2_0 v0 u(2)]';
    t_span = linspace(t0, u(1), 1000);
    [~, y] = ode45(@(t, y) dynamics(y, g), t_span, y0);

    ceq = y(end,2);
    dydt_final = dynamics(y(end,:)', g);

    pT = [0 -1 0 0];
    [~, p] = ode45(@(t, p) p_ode(t, p, y, g, t_span), flip(t_span), pT);
    q = p(end,:);

    Dceq = [dydt_final(2); -q(4)];
end

function dpdt = p_ode(t, p, y, g, t_span)
    v = interp1(t_span, y(:,3), t);
    gamma = interp1(t_span, y(:,4), t);
    J = [0 0 0 0;
         0 0 0 0;
         cos(gamma) sin(gamma) 0 (g/v^2)*cos(gamma);
         -v*sin(gamma) v*cos(gamma) -g*cos(gamma) (g/v)*sin(gamma)];
    dpdt = -J*p;
end

function dydt = dynamics(y, g)
    v = y(3);
    gamma = y(4);
    dydt = [
        v*cos(gamma);
        v*sin(gamma);
        -g*sin(gamma);
        -(g/v)*cos(gamma)
    ];
end
