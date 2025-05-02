function main_constraint_drag
    persistent u_history J_history;
    clc; close all

    g = 9.81; % Gravity [m/s^2]
    rho = 1.225; % Air density [kg/m^3]
    Sw = 2; % Wing/Aerodynamic surface [m^2]
    Cd0 = 0.01; % Drag coefficient
    m = 1; % Object mass [kg]
    k = 0.05;
    Clalpha = 2*pi;
    Cl0 = 0;

    N = 500; % Discretization

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

    Cl = @(alpha) Cl0 + Clalpha*alpha;
    Cd = @(alpha) Cd0 + k*Cl(alpha).^2;
    D = @(V,alpha) 0.5*rho*Sw*V.^2.*Cd(alpha)/m;
   
    dCl = Clalpha;
    dCd = @(alpha) 2*k*Cl(alpha)*dCl;
    dDda = @(V,alpha) 0.5*rho*Sw*V.^2.*dCd(alpha)/m;

    dDragdv = @(V,alpha) rho*Sw*V.*Cd(alpha)/m;


    cost = @(u) f_cost(u, g, t0, v0, gamma0, x1_0, x2_0, D,dDda,dDragdv, N);
    constraint = @(u) groundConstraint(u, g, t0, v0, x1_0, x2_0, gamma0, D,dDda,dDragdv, N);

    [valid,err] = checkGradients(cost, u0, "Display","on")
    [valid2,err2] = checkGradients(constraint, u0, "Display","on")

    options = optimoptions("fmincon", ...
        "OutputFcn", @store_fmincon, ...
        "Algorithm", "sqp", ... 
        "DerivativeCheck", "on","Display","iter");

    [u_opt, ~] = fmincon(cost, u0, [], [], [], [], lb, ub, @(u) nonlcon(u,constraint), options);

    y0 = [x1_0 x2_0 v0 gamma0];
    t_span = linspace(t0, u_opt(1), N);
    [~, y] = ode45(@(t, y) dynamics(t, y, tf, u_opt(2:N+1), g, D, N), t_span, y0);

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

function [J, gradJ] = f_cost(u, g, t0, v0, gamma0, x1_0, x2_0,Drag,dDragdv,dDragda, N)
    y0 = [x1_0 x2_0 v0 gamma0];
    tf = u(1);
    alpha = u(2:N+1);
    t_span = linspace(t0, u(1), N);    
    [t, y] = ode45(@(t, y) dynamics(t, y, tf, alpha, g, Drag, N), t_span, y0);

    J = -y(end,1);
    dydt_final = dynamics(t_span(end), y(end,:)', tf, alpha, g, Drag, N);

    pT = [1 0 0 0];
    [~, p] = ode45(@(t, p) adjoint(t, p, y, g, t_span, alpha, dDragdv), flip(t_span), pT);
    

    ind = round((t - t_span(1)) / (t_span(2) - t_span(1))) + 1;
    v = y(ind,3);    
    %v = interp1(t_span, y(:,3), t);
    DFdu = -dDragda(v,alpha');

    gradJ = [-dydt_final(1); -DFdu.*p(:,3)];
end

function [ceq, Dceq] = groundConstraint(u, g, t0, v0, x1_0, x2_0, gamma0, Drag,dDragdv,dDragda, N)
    y0 = [x1_0 x2_0 v0 gamma0]';
    tf = u(1);
    t_span = linspace(t0, tf, N);
    alpha = u(2:N+1);
    [t, y] = ode45(@(t, y) dynamics(t, y, tf, alpha, g, Drag, N), t_span, y0);

    ceq = y(end,2);
    dydt_final = dynamics(t_span(end), y(end,:)', tf, alpha, g, Drag, N);

    pT = [0 -1 0 0];
    [~, p] = ode45(@(t, p) adjoint(t, p, y, g, t_span, alpha, dDragdv), flip(t_span), pT);
   

    %v = interp1(t_span, y(:,3), t);
    ind = round((t - t_span(1)) / (t_span(2) - t_span(1))) + 1;
    v = y(ind,3);    
    DFdu = dDragda(v,alpha');

    Dceq = [dydt_final(2); -DFdu.*p(:,3)];
end

function dpdt = adjoint(t, p, y, g, t_span, alpha, dDragdv)
    %v = interp1(t_span, y(:,3), t);
    %gamma = interp1(t_span, y(:,4), t);
    %alpha = interp1(t_span, alpha, t);

    ind = round((t - t_span(1)) / (t_span(2) - t_span(1))) + 1;
    v = y(ind,3);
    gamma = y(ind,4);    
    alpha = alpha(ind);

    J = [0 0 0 0;
         0 0 0 0;
         cos(gamma) sin(gamma) -dDragdv(v,alpha') (g/v^2)*cos(gamma);
         -v*sin(gamma) v*cos(gamma) -g*cos(gamma) (g/v)*sin(gamma)];
    dpdt = -J*p;
end

function dydt = dynamics(t, y, tf, alpha, g, Drag, N)
    v = y(3);
    gamma = y(4);

    t_span = linspace(0, tf, N);    
    ind = round((t - t_span(1)) / (t_span(2) - t_span(1))) + 1;
    alpha = alpha(ind);

    
    %alpha = interp1(t_span, alpha, t);

    dydt = [
        v*cos(gamma);
        v*sin(gamma);
        -g*sin(gamma)-Drag(v,alpha);
        -(g/v)*cos(gamma)
    ];
end
