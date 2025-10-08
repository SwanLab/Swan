close all;
clear;
clc;

% Load trained network
load('Tutorials/ChomogNetworkTutorial/Networks/network_stiffness_matrix.mat')

%% Initialize the optimization problem

% Define the problem to solve
studyCases = {'maxStiffness';
              'minStiffness';
              'maxConditionNumber';
              'minConditionNumber';
              'targetStiffness';
              'maxDiagonal';
              'minDiagonal'};

studyType = char(studyCases(1));  % maxStiffness

% Define the problem constraint parameters
lower_bounds = 0.1;
upper_bounds = 2.0;

% Set the initial guess
x0 = 1.0;

% Wrap the cost and optimization functions
costfunction = @(x) stiffnessMatrixCost(opt, studyType, x);
nonlinconstraint = @(x) stiffnessMatrixConstraint(x);

% Set optimization problem options
options = optimoptions('fmincon', ...
                       'StepTolerance', 1.0e-16, ...
                       'SpecifyObjectiveGradient',true, ...
                       'SpecifyConstraintGradient',true, ...
                       'Display', 'iter', ...
                       'Algorithm', 'interior-point');

% Call optimizer
[x_opt, fval] = fmincon(costfunction, x0, [], [], [], [], lower_bounds, upper_bounds, nonlinconstraint, options);

%% Display results

% Get the optimized stiffness matrix
K_opt_vector = opt.computeOutputValues(x_opt);
K_opt = StiffnessMatrixHandler.vectorToMatrix(K_opt_vector);

disp('Optimal input parameter:');
disp(x_opt);
disp('Optimized cost:');
disp(fval);

disp('Optimized Stiffness Matrix (8x8):')
disp(K_opt)

% Analyze matrix properties
[condition_number, max_eigenval, min_eigenval] = StiffnessMatrixHandler.analyzeMatrix(K_opt);
fprintf('Matrix Properties:\n');
fprintf('  Condition Number: %.2e\n', condition_number);
fprintf('  Max Eigenvalue: %.2e\n', max_eigenval);
fprintf('  Min Eigenvalue: %.2e\n', min_eigenval);
fprintf('  Is Symmetric: %s\n', mat2str(issymmetric(K_opt, 1e-10)));
fprintf('  Is Positive Definite: %s\n', mat2str(all(eig(K_opt) > 0)));

% Visualize the optimized matrix
StiffnessMatrixHandler.plotStiffnessMatrix(K_opt, 'Optimized Stiffness Matrix');

%% Compare with different input parameters
figure;
tiledlayout(2, 2);

input_test = [0.5, 1.0, 1.5, 2.0];
for i = 1:4
    nexttile;
    
    K_test_vector = opt.computeOutputValues(input_test(i));
    K_test = StiffnessMatrixHandler.vectorToMatrix(K_test_vector);
    
    imagesc(K_test);
    colorbar;
    title(sprintf('Input = %.1f', input_test(i)));
    xlabel('DOF');
    ylabel('DOF');
end

sgtitle('Stiffness Matrices for Different Input Parameters');

%% Declare optimization functions

function [J, grad] = stiffnessMatrixCost(opt_stiffness, type, x)
    % Calcula el costo y gradiente para optimización de matriz de rigidez
    % INPUT: opt_stiffness - objeto de red neuronal entrenada
    %        type - tipo de optimización
    %        x - parámetro de entrada
    % OUTPUT: J - valor del costo
    %         grad - gradiente del costo

    % Fetch output and jacobian of the network
    Y = opt_stiffness.computeOutputValues(x);
    dY = opt_stiffness.computeGradient(x);
    
    % Reshape to 8x8 matrix
    K = StiffnessMatrixHandler.vectorToMatrix(Y);
    dK = reshape(dY, 8, 8);
    
    % Calculate the cost and gradient based on type
    switch type
        case 'maxStiffness'
            % Maximizar la traza (suma de elementos diagonales)
            J = -trace(K);
            grad = -trace(dK);
            
        case 'minStiffness'
            % Minimizar la traza
            J = trace(K);
            grad = trace(dK);
            
        case 'maxConditionNumber'
            % Maximizar el número de condición
            [condition_number, ~, ~] = StiffnessMatrixHandler.analyzeMatrix(K);
            J = -condition_number;
            % Gradiente aproximado (puede necesitar refinamiento)
            grad = -trace(dK) / trace(K) + condition_number * trace(dK) / trace(K);
            
        case 'minConditionNumber'
            % Minimizar el número de condición
            [condition_number, ~, ~] = StiffnessMatrixHandler.analyzeMatrix(K);
            J = condition_number;
            grad = trace(dK) / trace(K) - condition_number * trace(dK) / trace(K);
            
        case 'targetStiffness'
            % Acercarse a una matriz objetivo
            K_target = eye(8);  % Matriz objetivo (identidad)
            J = norm(K - K_target, 'fro')^2;
            grad = 2 * trace((K - K_target)' * dK);
            
        case 'maxDiagonal'
            % Maximizar la suma de elementos diagonales
            J = -sum(diag(K));
            grad = -sum(diag(dK));
            
        case 'minDiagonal'
            % Minimizar la suma de elementos diagonales
            J = sum(diag(K));
            grad = sum(diag(dK));
            
        otherwise
            error('Tipo de optimización no reconocido: %s', type);
    end
end

function [c, ceq, gradc, gradceq] = stiffnessMatrixConstraint(x)
    % Define constraints for the optimization problem
    % INPUT: x - parámetro de entrada
    % OUTPUT: c - desigualdades (c <= 0)
    %         ceq - igualdades (ceq = 0)
    %         gradc - gradiente de desigualdades
    %         gradceq - gradiente de igualdades

    % No inequality constraints
    c = [];
    gradc = [];
    
    % No equality constraints
    ceq = [];
    gradceq = [];
    
    % Ejemplo de restricción: asegurar que el parámetro esté en un rango específico
    % (esto ya está manejado por los bounds en fmincon)
    
    % Ejemplo de restricción adicional: asegurar que la matriz resultante sea válida
    % (esto se puede implementar si es necesario)
end
