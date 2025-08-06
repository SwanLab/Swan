close all;
clear;
clc;

% Load trained network
load('Tutorials/ChomogNetworkTutorial/Networks/network_ellipse_cHomog.mat')

%% Initialize the optimization problem

% Define the problem to solve
studyCases = {'maxHorzStiffness';
              'maxIsoStiffness';
              'maxAuxetic'};

studyType = char(studyCases(1));

% Define the problem constraint parameters
A_ellipse = 0.15^2 * pi;
lower_bounds = [10^-8; 10^-8];
upper_bounds = [0.5; 0.5];

% Set the initial guess
x0 = [rand(1) / 2; rand(1) / 2];

% Wrap the cost and optimization functions
costfunction = @(x) cHomogCost(opt, studyType, x);
nonlinconstraint = @(x) volumeConstraint(A_ellipse, x);

% Set optimization problem options
options = optimoptions('fmincon', 'StepTolerance', 1.0e-14,'SpecifyObjectiveGradient',true, 'SpecifyConstraintGradient',true, 'Display', 'iter');

% Call optimizer
[x_opt, fval] = fmincon(costfunction, x0, [], [], [], [], lower_bounds, upper_bounds, nonlinconstraint, options);

%% Display results

disp('Optimal solution:');
disp(x_opt);
disp('Optimal cost:');
disp(fval);
drawRectWithEllipseHole([-0.5, -0.5, 1, 1], [0, 0, x_opt(1), x_opt(2)])

Chomog = opt.computeOutputValues(x_opt');
poisson = Chomog(2) / Chomog(1);
disp('Poisson ratio:')
disp(poisson)

%% Declare optimization functions

function [c, ceq, gradc, gradceq] = volumeConstraint(A_ellipse, x)

    % No inequality constraints
    c = [];

    % Equality constraint: x1 - k/x2 == 0
    k = A_ellipse / pi;
    ceq = x(1) - k / x(2);

    if nargout > 2
        % Gradient of inequality constraints (none here)
        gradc = [];

        % Gradient of equality constraint
        % d/dx1 = 1
        % d/dx2 = k / x2^2
        gradceq = [1; k / (x(2)^2)];
    end
end

function drawRectWithEllipseHole(rect, ellipse)
% Draws a black rectangle with a white elliptical hole
% 
% rect = [x, y, width, height]    — rectangle lower-left corner and size
% ellipse = [cx, cy, rx, ry]      — ellipse center and radii

    % Unpack rectangle
    x = rect(1);
    y = rect(2);
    w = rect(3);
    h = rect(4);
    
    % Unpack ellipse
    cx = ellipse(1);
    cy = ellipse(2);
    rx = ellipse(3);
    ry = ellipse(4);
    
    % Draw black rectangle
    fill([x, x+w, x+w, x], [y, y, y+h, y+h], 'k');
    hold on;

    % Create ellipse coordinates
    theta = linspace(0, 2*pi, 200);
    ex = cx + rx * cos(theta);
    ey = cy + ry * sin(theta);

    % Draw white ellipse over the rectangle (to make a "hole")
    fill(ex, ey, 'w');

    axis equal;
    xlim([x, x+w]);
    ylim([y, y+h]);
    set(gca, 'Color', 'w');

    xlabel('x')
    ylabel('y')
    title(['Ellipse area: ', num2str(pi * rx * ry / (h * w) * 100, '%.1f'), ' % of rectangle area'])
    hold off;
end
