close all;
clear;
clc;

% Load trained networks
load('Tutorials/ChomogNetworkTutorial/Networks/network_superForm_area.mat')
load('Tutorials/ChomogNetworkTutorial/Networks/network_superForm_cHomog.mat')
sf = superformula_functionality;

%% Initialize the optimization problem

% Define the problem to solve
studyCases = {'maxHorzStiffness';
              'minHorzStiffness';
              'maxIsoStiffness';
              'maxAuxetic'};

studyType = char(studyCases(2));

% Define the problem constraint parameters - X vector: a, n (n1 = n2 = n3)
A_target = 0.2;
lower_bounds = [0.2; 2];
upper_bounds = [0.4; 12];

% Set the initial guess
x0 = [lower_bounds(1) + rand(1) * (upper_bounds(1) - lower_bounds(1));
      lower_bounds(2) + rand(1) * (upper_bounds(2) - lower_bounds(2))];

% Wrap the cost and optimization functions
costfunction = @(x) cHomogCost(opt_cHomog, studyType, x);
nonlinconstraint = @(x) volumeConstraint(A_target, sf, opt_area, x);

% Set optimization problem options
options = optimoptions('fmincon', 'StepTolerance', 1.0e-14,'SpecifyObjectiveGradient',true, 'SpecifyConstraintGradient',true, 'Display', 'iter');

% Call optimizer
[x_opt, fval] = fmincon(costfunction, x0, [], [], [], [], lower_bounds, upper_bounds, nonlinconstraint, options);
%[x_opt, fval] = fmincon(costfunction, x0, [], [], [], [], lower_bounds, upper_bounds, [], options);

%% Display results

gPar.semiVerticalAxis = x_opt(1); %b
gPar.m  = 8;
gPar.n1 = x_opt(2);
gPar.n2 = x_opt(2);
gPar.n3 = x_opt(2);
gPar.semiHorizontalAxis = gPar.semiVerticalAxis^(gPar.n3/gPar.n2); % S'ajusta a per un b donat per aconseguir un radi desitjat!!

rectangle = [-0.5, -0.5, 1, 1];

drawRectWithSuperHole(rectangle, gPar, sf)

disp('Optimal solution:');
disp(x_opt);
disp('Optimal cost:');
disp(fval);

Chomog = opt_cHomog.computeOutputValues(x_opt');
poisson = Chomog(2) / Chomog(1);
disp('Poisson ratio:')
disp(poisson)

%% Declare optimization functions

function [J, grad] = cHomogCost(opt_cHomog, type, x)

    % Fetch output and jacobian of the network
    Y = opt_cHomog.computeOutputValues(x');
    dY = opt_cHomog.computeGradient(x');
    
    % Fetch the Chomog component to optimize
    C_11 = Y(1);
    C_12 = Y(2);
    C_22 = Y(3);
    dC_11 = dY(:, 1)';
    dC_12 = dY(:, 2)';
    dC_22 = dY(:, 3)';

    % Calculate the cost and gradient
    switch type
        case 'maxHorzStiffness'
            J = 1 / C_11;
            grad = - 1 / C_11^2 .* dC_11;
        case 'minHorzStiffness'
            J = C_11;
            grad = dC_11;
        case 'maxIsoStiffness'
            J = 1 / (C_11 + C_22);
            grad = - (dC_11 + dC_22) / (C_11 + C_22)^2;
        case 'maxAuxetic'
            J = C_12 / C_11;
            grad = dC_12 / C_11 - C_12 / C_11^2 * dC_11;
    end

end

function [c, ceq, gradc, gradceq] = volumeConstraint(A_target, sf, A_net, x)

    % No inequality constraints
    c = [];

    % Equality constraint: volume residue
    b = x(1); %b
    m  = 8;
    n1 = x(2);
    n2 = x(2);
    n3 = x(2);
    a = b^(n3/n2); % S'ajusta a per un b donat per aconseguir un radi desitjat!!

    A_current = sf.compute_area(a, b, m, n1, n2, n3);
    dA_current = A_net.computeGradient([a, n1]);

    ceq = A_target - A_current;

    if nargout > 2
        % Gradient of inequality constraints (none here)
        gradc = [];

        % Gradient of equality constraint
        gradceq = - [dA_current(1); dA_current(2)];

    end
end

function drawRectWithSuperHole(rect, gPar, sf)

    x = rect(1);
    y = rect(2);
    w = rect(3);
    h = rect(4);

    a = gPar.semiHorizontalAxis;
    b = gPar.semiVerticalAxis;   
    
    m = gPar.m;
    n1 = gPar.n1;
    n2 = gPar.n2;
    n3 = gPar.n3;

    area_sup = sf.compute_area(a, b, m, n1, n2, n3);

    phi_vec = linspace(0, 2*pi, 1000);
    rad_vec = sf.calculate(phi_vec, a, b, m, n1, n2, n3);
    x_vec = rad_vec .* cos(phi_vec);
    y_vec = rad_vec .* sin(phi_vec);

    % Draw black rectangle
    fill([x, x+w, x+w, x], [y, y, y+h, y+h], 'k');
    hold on;

    % Draw white ellipse over the rectangle (to make a "hole")
    fill(x_vec, y_vec, 'w');

    axis equal;
    xlim([x, x+w]);
    ylim([y, y+h]);
    set(gca, 'Color', 'w');

    xlabel('x')
    ylabel('y')
    title(['Ellipse area: ', num2str(area_sup / (w * h) * 100 , '%.1f'), ' % of rectangle area'])
    hold off;

    axis equal

end