%% run_Al7075_WF110_cases.m — Batch script for Al7075_WF110
%  Runs STATIC and DYNAMIC cases, saves figures and console output.

addpath(genpath('.'));

resultsRoot = 'RESULTS_Al7075_WF110';
if ~exist(resultsRoot, 'dir')
    mkdir(resultsRoot);
end

totalCases = 0;
skippedCases = 0;
failedCases = 0;
casesRun = 0;

fprintf('==========================================================\n');
fprintf('  AUTOMATIC OPTIMIZATION RUNNER — Started %s\n', datestr(now));
fprintf('==========================================================\n\n');

%% ========================================================================
%  Case 1: STATIC, forceSymmetry = true
%% ========================================================================
totalCases = totalCases + 1;
caseName = 'STATIC_sym_on';
caseDir = fullfile(resultsRoot, caseName);
if exist(fullfile(caseDir, 'output.txt'), 'file')
    fprintf('[SKIP] %s (already exists)\n', caseName);
    skippedCases = skippedCases + 1;
else
    fprintf('[RUN]  %s ...\n', caseName);
    if ~exist(caseDir, 'dir'), mkdir(caseDir); end
    diary(fullfile(caseDir, 'output.txt'));
    try
        Al7075_WF110_('STATIC', true);
        saveFigures(caseDir);
        casesRun = casesRun + 1;
    catch ME
        fprintf('ERROR: %s\n', ME.message);
        disp(getReport(ME));
        failedCases = failedCases + 1;
    end
    close all; diary off;
end

%% ========================================================================
%  Case 2: DYNAMIC, forceSymmetry = true
%% ========================================================================
totalCases = totalCases + 1;
caseName = 'DYNAMIC_sym_on';
caseDir = fullfile(resultsRoot, caseName);
if exist(fullfile(caseDir, 'output.txt'), 'file')
    fprintf('[SKIP] %s (already exists)\n', caseName);
    skippedCases = skippedCases + 1;
else
    fprintf('[RUN]  %s ...\n', caseName);
    if ~exist(caseDir, 'dir'), mkdir(caseDir); end
    diary(fullfile(caseDir, 'output.txt'));
    try
        Al7075_WF110_('DYNAMIC', true);
        saveFigures(caseDir);
        casesRun = casesRun + 1;
    catch ME
        fprintf('ERROR: %s\n', ME.message);
        disp(getReport(ME));
        failedCases = failedCases + 1;
    end
    close all; diary off;
end

%% ========================================================================
%  SUMMARY
%% ========================================================================
fprintf('\n==========================================================\n');
fprintf('  RUNNER SUMMARY — Finished %s\n', datestr(now));
fprintf('==========================================================\n');
fprintf('  Total cases   : %d\n', totalCases);
fprintf('  Run           : %d\n', casesRun);
fprintf('  Skipped       : %d\n', skippedCases);
fprintf('  Failed        : %d\n', failedCases);
fprintf('==========================================================\n\n');

%% ========================================================================
%  HELPER FUNCTION — save all open figures as .fig and .png
%% ========================================================================
function saveFigures(caseDir)
    figs = findobj('Type', 'figure');
    for i = 1:length(figs)
        figName = sprintf('fig%d', figs(i).Number);
        saveas(figs(i), fullfile(caseDir, [figName '.png']));
        saveas(figs(i), fullfile(caseDir, [figName '.fig']));
    end
end
