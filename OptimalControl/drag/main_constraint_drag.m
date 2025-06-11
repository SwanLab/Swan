function main_constraint_drag
    persistent u_history J_history;
    clc; close all

    g = 9.81; % Gravity [m/s^2]
    rho = 1.225; % Air density [kg/m^3]
    Sw = 0.6; % Wing/Aerodynamic surface [m^2]
    Cd0 = 0.03; % Drag coefficient
    m = 5; % Object mass [kg]
    k = 0.05;
    Clalpha = 5.7296;
    Cl0 = 0;

    N = 50; % Discretization

    x1_0 = 0; x2_0 = 0; v0 = 15;
    gamma0 = deg2rad(40);
    t0 = 0; tf = 3;
    alpha0 = deg2rad(3);
    u0 = [tf; alpha0];
    lb = [0 0]; % Lower bounds for the control
    ub = [5 deg2rad(10)]; % Upper bounds for the control

    u0 = [u0(1) ones(1,N)*u0(2)];
    lb = [lb(1) ones(1,N)*lb(2)];
    ub = [ub(1) ones(1,N)*ub(2)];

    % Previous calculations (Derivatives)

    Cl = @(alpha) Cl0 + Clalpha*alpha;
    Cd = @(alpha) Cd0 + k*Cl(alpha).^2;
    D = @(V,alpha) 0.5*rho*Sw*V.^2.*Cd(alpha);
   
    dCl = Clalpha; % Derivative respecto to alpha
    dCd = @(alpha) 2*k*Cl(alpha)*dCl; % Derivative respect to alpha
    dDda = @(V,alpha) 0.5*rho*Sw*V.^2.*dCd(alpha); % Derivative respect to alpha

    dDdv = @(V,alpha) rho*Sw*V.*Cd(alpha); % Derivative respect to velocity

    cost = @(u) f_cost(u, g, t0, v0, gamma0, x1_0, x2_0, D, dDda, dDdv, N, m);
    constraint = @(u) groundConstraint(u, g, t0, v0, x1_0, x2_0, gamma0, D, dDda, dDdv, N, m);
    % 
    % eps = 1e-6;
    % [cost0, grad0] = cost(u0);
    % for i=1:N+1
    %     uI = u0;
    %     uI(i) = uI(i) + eps*u0(i);
    %     costI(i) = cost(uI);
    %     gradD(i) = (costI(i)-cost0)/(uI(i)-u0(i));
    % end

    % plot(grad0(1:end),'-+')
    % hold on
    % plot(gradD(1:end),'-+')

    % Coste centrada

    % eps = 1e-6;
    % [cost0, grad0] = cost(u0);
    % 
    % gradD = zeros(N+1, 1);  % inicializar vector de derivadas
    % 
    % for i = 1:N+1
    %     u_forward = u0;
    %     u_backward = u0;
    % 
    %     % Asegura que el paso no sea nulo si u0(i) es cero
    %     delta = eps * max(1, abs(u0(i)));
    % 
    %     u_forward(i) = u_forward(i) + delta*u0(i);
    %     u_backward(i) = u_backward(i) - delta*u0(i);
    % 
    %     cost_forward = cost(u_forward);
    %     cost_backward = cost(u_backward);
    % 
    %     % Diferencia centrada
    %     gradD(i) = (cost_forward - cost_backward) / (2 * delta);
    % end
    % 
    % plot(grad0(2:end)/(tf/N), '-+')
    % hold on
    % plot(gradD(2:end), '-+')
    % legend('Gradiente analítico', 'Gradiente centrado numérico')
    % xlabel('Índice de variable de control')
    % ylabel('Derivada')
    % title('Comparación entre gradiente del coste analítico y numérico centrado')


    % eps = 1e-6;
    % [constraint0,gradconst0] = constraint(u0);
    % for i=1:N+1
    %     uI = u0;
    %     uI(i) = uI(i) + eps*u0(i);
    %     constraintI(i) = constraint(uI);
    %     gradconstD(i) = (constraintI(i)-constraint0)/(uI(i)-u0(i));
    % end
    % 
    % figure
    % plot(gradconst0/(tf/N),'-+')
    % hold on
    % plot(gradconstD,'-+')

    % Centrada

    % eps = 1e-6;
    % [constraint0, gradconst0] = constraint(u0);
    % 
    % gradconstD = zeros(N+1, 1);  % inicializar vector de derivadas
    % 
    % for i = 1:N+1
    %     u_forward = u0;
    %     u_backward = u0;
    % 
    %     % evitar que eps*u0(i) sea cero si u0(i) lo es
    %     delta = eps * max(1, abs(u0(i)));
    % 
    %     u_forward(i) = u_forward(i) + delta;
    %     u_backward(i) = u_backward(i) - delta;
    % 
    %     constraint_forward = constraint(u_forward);
    %     constraint_backward = constraint(u_backward);
    % 
    %     % derivada centrada
    %     gradconstD(i) = (constraint_forward - constraint_backward) / (2 * delta);
    % end
    % 
    % figure
    % plot(gradconst0(2:end)/(tf/N), '-+')
    % hold on
    % plot(gradconstD(2:end), '-+')
    % legend('Gradiente analítico', 'Gradiente centrado numérico')
    % xlabel('Índice de variable de control')
    % ylabel('Derivada')
    % title('Comparación entre gradiente analítico y numérico centrado')
    % 
    % [valid,err] = checkGradients(cost, u0, 'Tolerance', 1e-3, 'Display', 'on');
    % assignin('base', 'err', err);
    % [valid2,err2] = checkGradients(constraint, u0, 'Tolerance', 1e-3, 'Display', 'on')
    % assignin('base', 'err2', err2);


    options = optimoptions("fmincon", ...
        "OutputFcn", @store_fmincon, ...
        "Algorithm", "sqp", ... 
        "DerivativeCheck", "on","Display","iter");

    [u_opt, ~] = fmincon(cost, u0, [], [], [], [], lb, ub, @(u) nonlcon(u,constraint), options);

    y0 = [x1_0 x2_0 v0 gamma0];
    t_span = linspace(t0, u_opt(1), N);
    [~, y] = ode45(@(t, y) dynamics(t, y, tf, u_opt(2:N+1), g, D, N, m), t_span, y0);

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

function [c, ceq] = nonlcon(u,constraint)
    c = [];
    Dc = [];
    [ceq] = constraint(u);
end

function [J] = f_cost(u, g, t0, v0, gamma0, x1_0, x2_0, D, dDdv, dDda, N, m)
    y0 = [x1_0 x2_0 v0 gamma0];
    tf = u(1);
    alpha = u(2:N+1);
    t_span = linspace(t0, tf, N)';    
    [~, y] = ode45(@(t, y) dynamics(t, y, tf, alpha, g, D, N, m), t_span, y0);

    J = -y(end,1);
    dydt_final = dynamics(tf, y(end,:)', tf, alpha, g, D, N, m);

    pT = [1 0 0 0];
    [~, p] = ode45(@(t, p) adjoint(t, p, y, g, t_span, alpha, dDdv, m), flip(t_span), pT);

    v = y(:,3);
    DFdu = -dDda(v,alpha')./m;
    p3 = p((N:-1:1),3);
    gradJ = [-dydt_final(1); -DFdu.*p3];


    figure(100); plot(y(:,1), y(:,2), 'b-', 'LineWidth', 2);
    xlabel("Horizontal distance [m]"); ylabel("Vertical distance [m]"); grid on;
    xlim([0 30])
    ylim([-1 5])

end

function [ceq] = groundConstraint(u, g, t0, v0, x1_0, x2_0, gamma0, D, dDdv, dDda, N, m)
    y0 = [x1_0 x2_0 v0 gamma0]';
    tf = u(1);
    t_span = linspace(t0, tf, N);
    alpha = u(2:N+1);
    [~, y] = ode45(@(t, y) dynamics(t, y, tf, alpha, g, D, N, m), t_span, y0);

    ceq = y(end,2);
    dydt_final = dynamics(t_span(end), y(end,:)', tf, alpha, g, D, N, m);

    pT = [0 -1 0 0];
    [~, p] = ode45(@(t, p) adjoint(t, p, y, g, t_span, alpha, dDdv, m), flip(t_span), pT);
     
    v = y(:,3);
    DFdu = -dDda(v,alpha');

    Dceq = [dydt_final(2); -DFdu.*p(:,3)];
end

function dpdt = adjoint(t, p, y, g, t_span, alpha, dDdv, m)
    v = interp1(t_span, y(:,3), t);
    gamma = interp1(t_span, y(:,4), t);
    alpha = interp1(t_span, alpha, t);

    J = [0 0 0 0;
         0 0 0 0;
         cos(gamma) sin(gamma) -dDdv(v,alpha)./m (g/v^2)*cos(gamma);
         -v*sin(gamma) v*cos(gamma) -g*cos(gamma) (g/v)*sin(gamma)];
    dpdt = -J*p;
end

function dydt = dynamics(t, y, tf, alpha, g, D, N, m)
    v = y(3);
    gamma = y(4);

    t_span = linspace(0, tf, N); 
    alpha = interp1(t_span, alpha, t);

    dydt = [
        v*cos(gamma);
        v*sin(gamma);
        -g*sin(gamma)-D(v,alpha)/m;
        -(g/v)*cos(gamma)
    ];
end
