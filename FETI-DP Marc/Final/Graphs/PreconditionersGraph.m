%% Script for FETI-DP Preconditioners Convergence Comparison
clear; clc; close all;

fprintf('===================================================\n');
fprintf('     STARTING PRECONDITIONER COMPARISON ANALYSIS  \n');
fprintf('===================================================\n');

% Preconditioners configuration
preconditioners = {'dirichlet', 'ilu', 'amg'};
colors = { ...
    [0.47, 0.67, 0.19], ... % Green for Dirichlet
    [0.93, 0.69, 0.13], ... % Yellow/Amber for ILU
    [0.64, 0.08, 0.18]      % Crimson for AMG
};
markers = {'-^', '-v', '-s'};

% Structures to store iteration histories
fetiResiduals = struct();
legendLabels = {};

% Flag to compute and extract the unpreconditioned system history once
saveUnpreconditioned = true;
unprecondResidual = [];

%% 1. Case Execution Loop
for i = 1:length(preconditioners)
    currentPrec = preconditioners{i};
    
    fprintf('\n>>> Executing case with preconditioner: %s\n', upper(currentPrec));
    
    % Instantiate class using the updated constructor
    runner = DataGenerationElasticity2DPreconditioners(currentPrec);
    
    % Override properties to ensure a clean benchmark run
    runner.enablePlots             = false; % Suppress individual system plots
    runner.computeMonolithic       = false; % Disable to accelerate the loop
    runner.computeUnpreconditioned = saveUnpreconditioned; 
    
    % Compute dual interface problem and extract convergence history
    results = runner.solveFetiDP();
    
    % Store the preconditioned history
    fetiResiduals.(currentPrec) = results.residual;
    
    % Capture unpreconditioned data from the first run only
    if saveUnpreconditioned && isfield(results, 'residualUnprec')
        unprecondResidual = results.residualUnprec;
        saveUnpreconditioned = false; 
    end
end

%% 2. Generation of the Final Publication-Ready Figure
figure('Name', 'FETI-DP Preconditioner Convergence Comparison', ...
       'Color', 'w', 'Position', [100, 100, 850, 520]);
hold on;

% 2.1. Plot Unpreconditioned Case (FETI dual CG)
if ~isempty(unprecondResidual)
    semilogy(1:length(unprecondResidual), unprecondResidual, '-o', ...
        'Color', [0.85, 0.33, 0.10], 'LineWidth', 2, 'MarkerSize', 4);
    legendLabels{end+1} = 'FETI-DP CG';
end

% 2.2. Plot Preconditioned Cases
for i = 1:length(preconditioners)
    currentPrec = preconditioners{i};
    res = fetiResiduals.(currentPrec);
    
    semilogy(1:length(res), res, markers{i}, ...
        'Color', colors{i}, 'LineWidth', 2, 'MarkerSize', 5);
    legendLabels{end+1} = sprintf('FETI-DP PCG (%s)', upper(currentPrec));
end

% 2.3. Convergence Tolerance Line
tolVal = runner.pcgTol;
yline(tolVal, '--k', 'LineWidth', 1.5, ...
      'Label', sprintf('Tolerance = %.0e', tolVal), ...
      'LabelHorizontalAlignment', 'left', 'FontSize', 10);

%% 3. Figure Formatting and Aesthetics
grid on;
ax = gca;
ax.GridAlpha = 0.25;
ax.YMinorGrid = 'on';
ax.YScale = 'log';

% Academic formatting labels
xlabel('Iteration Count', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Relative Residual \mid\midr_k\mid\mid_2 / \mid\midr_0\mid\mid_2', 'FontSize', 12, 'FontWeight', 'bold');
title('Dual Conjugate Gradient Convergence History - 2D Elasticity', 'FontSize', 13, 'FontWeight', 'bold');

legend(legendLabels, 'Location', 'northeast', 'FontSize', 11, 'Box', 'on');

% Set tight, dynamic axis boundaries
maxIter = length(unprecondResidual);
for i = 1:length(preconditioners)
    maxIter = max(maxIter, length(fetiResiduals.(preconditioners{i})));
end
xlim([1, maxIter + 2]);
ylim([tolVal * 0.1, 2]);

hold off;
fprintf('\n===================================================\n');
fprintf('     COMPARISON GRAPH GENERATED SUCCESSFULY\n');
fprintf('===================================================\n');