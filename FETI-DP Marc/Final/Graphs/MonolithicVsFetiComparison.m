%% Script for Comprehensive Convergence Comparison (Monolithic vs. FETI-DP)
clear; clc; close all;

fprintf('===================================================\n');
fprintf('   STARTING MONOLITHIC VS. FETI-DP BENCHMARK      \n');
fprintf('===================================================\n');

% Monolithic configuration cases to loop through (Lengths matched to 2)
monoPreconditioners = {'ilu', 'amg'}; 
monoColors = { ...
    [0.93, 0.69, 0.13], ... % Amber for Monolithic PCG (ILU)
    [0.64, 0.08, 0.18]      % Crimson for Monolithic PCG (AMG)
};
monoLabels = { ...
    'Monolithic PCG (ILU)', ...
    'Monolithic PCG (AMG)' ...
};

% Data allocation structures
monolithicResiduals = struct();
fetiDirichletRes = [];
fetiUnprecRes    = [];
legendLabels     = {};

% Tracking flags to extract shared FETI cases only once
captureFetiBaselines = true;

%% 1. Case Execution and Convergence Data Gathering
for i = 1:length(monoPreconditioners)
    currentMonoPrec = monoPreconditioners{i};
    
    fprintf('\n>>> Benchmarking Monolithic Preconditioner: %s\n', upper(currentMonoPrec));
    
    % FIXED: Passing false here prevents the constructor from auto-solving,
    % avoiding duplicate solver logs and double execution time.
    runner = DataGenerationElasticity2DMonoPrecond(currentMonoPrec, false);
    
    % Override flags via public properties before computing to guarantee a clean runtime
    runner.enablePlots             = false; % Suppress separate figure windows
    runner.computeMonolithic       = true;  % Must be true to extract monolithic data
    runner.computeUnpreconditioned = captureFetiBaselines; 
    
    % Execute solving pipeline (Runs exactly once with your custom flags)
    results = runner.solveFetiDP();
    
    % Save Monolithic residuals history dynamically
    monolithicResiduals.(currentMonoPrec) = results.residualMono;
    
    % Save the FETI interface data from the first run only
    if captureFetiBaselines
        if isfield(results, 'residual')
            fetiDirichletRes = results.residual;
        end
        if isfield(results, 'residualUnprec')
            fetiUnprecRes = results.residualUnprec;
        end
        captureFetiBaselines = false; % Locked for the following steps
    end
end

%% 2. Generation of the Publication-Ready Comparison Figure
figure('Name', 'Global Solver Convergence History Comparison', ...
       'Color', 'w', 'Position', [100, 100, 900, 550]);
hold on;

% 2.1. Plot Domain Decomposition (FETI) Curves first
if ~isempty(fetiUnprecRes)
    semilogy(1:length(fetiUnprecRes), fetiUnprecRes, '-s', ...
        'Color', [0.85, 0.33, 0.10], 'LineWidth', 2, 'MarkerSize', 4);
    legendLabels{end+1} = 'FETI-DP Unpreconditioned';
end

if ~isempty(fetiDirichletRes)
    semilogy(1:length(fetiDirichletRes), fetiDirichletRes, '-^', ...
        'Color', [0.47, 0.67, 0.19], 'LineWidth', 2, 'MarkerSize', 5);
    legendLabels{end+1} = 'FETI-DP PCG (Dirichlet)';
end

% 2.2. Plot Monolithic Convergence History Curves
for i = 1:length(monoPreconditioners)
    currentMonoPrec = monoPreconditioners{i};
    resMono = monolithicResiduals.(currentMonoPrec);
    
    if ~isempty(resMono)
        semilogy(1:length(resMono), resMono, '-o', ...
            'Color', monoColors{i}, 'LineWidth', 2, 'MarkerSize', 4);
        legendLabels{end+1} = monoLabels{i};
    end
end

% 2.3. Convergence Limit Reference Boundary
tolVal = runner.pcgTol;
yline(tolVal, '--k', 'LineWidth', 1.5, ...
      'Label', sprintf('Tolerance = %.0e', tolVal), ...
      'LabelHorizontalAlignment', 'left', 'FontSize', 10);

%% 3. Academic Layout Formatting
grid on;
ax = gca;
ax.GridAlpha = 0.25;
ax.YMinorGrid = 'on';
ax.YScale = 'log';

% Axis Titles and Mathematical Notation
xlabel('Iteration Count', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Relative Residual \mid\midr_k\mid\mid_2 / \mid\midr_0\mid\mid_2', 'FontSize', 12, 'FontWeight', 'bold');
title('Global Convergence Performance - Monolithic vs. Domain Decomposition (FETI-DP)', ...
      'FontSize', 12, 'FontWeight', 'bold');

legend(legendLabels, 'Location', 'northeast', 'FontSize', 10, 'Box', 'on');

% Dynamic bounding box limits calibration
maxIterationsObserved = max([length(fetiUnprecRes), length(fetiDirichletRes)]);
for i = 1:length(monoPreconditioners)
    maxIterationsObserved = max(maxIterationsObserved, length(monolithicResiduals.(monoPreconditioners{i})));
end
xlim([1, maxIterationsObserved + 5]);
ylim([tolVal * 0.1, 2]);

hold off;
fprintf('\n===================================================\n');
fprintf('     COMPARATIVE GRAPH EXPORTED SUCCESSFULY\n');
fprintf('===================================================\n');