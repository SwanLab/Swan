function main_constraint_draglift
    persistent u_history J_history;
    clc; close all

    data = load('resultsfalcon.mat');
    resultsfalcon = data.results;
    clear data
    g = 9.81; % Gravity [m/s^2]
    rho = 1.225; % Air density [kg/m^3]
    Sw = 0.6; % Wing/Aerodynamic surface [m^2]
    Cd0 = 0.03; % Drag coefficient
    m = 5; % Object mass [kg]
    k = 0.05;
    Clalpha = 0.1*360/(2*pi); % en radianes
    %Clalpha = 0.1; % en grados
    Cl0 = 0;

    N = 250; % Discretization

    x1_0 = 0; x2_0 = 0; v0 = 15;
    gamma0 = deg2rad(30);
    t0 = 0; tf_guess = 6;
    alpha0 = deg2rad(10);
    u0 = [resultsfalcon.tfinal resultsfalcon.u.alpha];
    %u0 = [tf_guess alpha0]
    lb = [0.1 deg2rad(-10)]; % Lower bounds for the control
    ub = [15 deg2rad(10.1)]; % Upper bounds for the control

    %u0 = [u0(1) ones(1,N)*u0(2)];
    lb = [lb(1) ones(1,N)*lb(2)];
    ub = [ub(1) ones(1,N)*ub(2)];

    % Previous calculations (Derivatives)

    Cl = @(alpha) Cl0 + Clalpha*alpha;
    Cd = @(alpha) Cd0 + k*Cl(alpha).^2;
    D = @(V,alpha) 0.5*rho*Sw*V.^2.*Cd(alpha);
    L = @(V,alpha) 0.5*rho*Sw*V.^2.*Cl(alpha);   

    [~,y0traj] = ode45(@(t,y)dynamics(t,y,u0(1),u0(2:end),g,D,N,m,L), ...
                   linspace(t0,u0(1),N), ...
                   [x1_0; x2_0; v0; gamma0]);

minAlt = min(y0traj(:,2));
fprintf('Altura mínima con u0 = %.4f m\n',minAlt);

    cost = @(u) f_cost(u, g, t0, v0, gamma0, x1_0, x2_0, D, N, m, L);
    nonlcon_fun = @(u) nonlcon_full(u, g, t0, v0, x1_0, x2_0, gamma0, D, N, m, L);

    options = optimoptions("fmincon", ...
        "OutputFcn", @store_fmincon, ...
        "Algorithm", "sqp", ... 
        "DerivativeCheck", "off","Display","iter","FiniteDifferenceType","central", ...
        "FiniteDifferenceStepSize", 1e-4);

    [u_opt, ~] = fmincon(cost, u0, [], [], [], [], lb, ub, nonlcon_fun, options);

    y0 = [x1_0 x2_0 v0 gamma0];
    tf = u_opt(1);
    t_span = linspace(t0, tf, N);
    [~, y] = ode45(@(t, y) dynamics(t, y, tf, u_opt(2:N+1), g, D, N, m, L), t_span, y0);

    disp(["Maximum distance [m] = ", num2str(y(end, 1))])
    disp(["Final Time [s]: ", num2str(u_opt(1))])

    figure
    plot(y(:,1),y(:,2))
    xlabel("Horizontal displacement [m]")
    ylabel("Vertical displacement [m]")
    title("Trajectory")

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

function [J] = f_cost(u, g, t0, v0, gamma0, x1_0, x2_0, D, N, m, L)
    y0 = [x1_0 x2_0 v0 gamma0];
    tf = u(1);
    alpha = u(2:N+1);
    t_span = linspace(t0, tf, N)';
    ode_opts = odeset('RelTol',1e-5,'AbsTol',1e-7);
    [~, y] = ode45(@(t, y) dynamics(t, y, tf, alpha, g, D, N, m, L), t_span, y0, ode_opts);

    J = -y(end,1);

    % figure(100); plot(y(:,1), y(:,2), 'b-+', 'LineWidth', 2);
    % xlabel("Horizontal distance [m]"); ylabel("Vertical distance [m]"); grid on;
    % xlim([0 200])
    % ylim([-50 50])
    % 
    % figure(101); plot(t_span, rad2deg(u(2:end)), 'b-+', 'LineWidth', 2);
    % xlabel("Time [t]"); ylabel("angle of attack [deg]"); grid on;
    % xlim([0 tf])
    % ylim([0 10])
    % 
    % figure(102); plot(t_span, (y(:,3)), 'b-+', 'LineWidth', 2);
    % xlabel("Time [t]"); ylabel("Velocity [m/s]"); grid on;
    % xlim([0 tf])
    % ylim([-50 50])
    % 
    % figure(103); plot(J);
    % hold on

end

function [c, ceq] = nonlcon_full(u, g, t0, v0, gamma0, x1_0, x2_0, D, N, m, L)
    tf    = u(1);
    alpha = u(2:end);
    t_span = linspace(t0, tf, N);
    y0 = [x1_0; x2_0; v0; gamma0];
    ode_opts = odeset('RelTol',1e-5,'AbsTol',1e-7);
    [~, y] = ode45(@(t,y)dynamics(t,y,tf,alpha,g,D,N,m,L), t_span, y0, ode_opts);

    c = -y(:,2);
    ceq = y(end,2);

end

function dydt = dynamics(t, y, tf, alpha, g, D, N, m, L)
    v = y(3);
    gamma = y(4);

    t_span = linspace(0, tf, N); 
    alpha = interp1(t_span, alpha, t);

    dydt = [
        v*cos(gamma);
        v*sin(gamma);
        -g*sin(gamma)-D(v,alpha)/m;
        -(g/v)*cos(gamma)+L(v,alpha)/(m*v);
    ];
end
