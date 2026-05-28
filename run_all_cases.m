%% run_all_cases.m — Master runner for all optimization cases
%  Runs all cases from CODES_TO_RUN.txt, saves figures (.fig/.png) and
%  command window output in RESULTS_T300/ with subfolders per code.
%  Skips cases whose results already exist (checks for output.txt).
%
%  Usage: matlab -batch "run('run_all_cases.m')"

addpath(genpath('.'));

resultsRoot = 'RESULTS_T300';
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
%  1. TutorialShellsOptim — STATIC
%     Case A: forceSymmetry = false
%     Case B: forceSymmetry = true
%% ========================================================================
codeFolder = fullfile(resultsRoot, 'TutorialShellsOptim');
if ~exist(codeFolder, 'dir'), mkdir(codeFolder); end

% --- Case 1A: static_sym_off ---
totalCases = totalCases + 1;
caseDir = fullfile(codeFolder, 'static_sym_off');
if exist(fullfile(caseDir, 'output.txt'), 'file')
    fprintf('[SKIP] TutorialShellsOptim / static_sym_off (already exists)\n');
    skippedCases = skippedCases + 1;
else
    fprintf('[RUN]  TutorialShellsOptim / static_sym_off ...\n');
    if ~exist(caseDir, 'dir'), mkdir(caseDir); end
    diary(fullfile(caseDir, 'output.txt'));
    try
        TutorialShellsOptim_(false);
        saveFigures(caseDir);
        casesRun = casesRun + 1;
    catch ME
        fprintf('ERROR: %s\n', ME.message);
        disp(getReport(ME));
        failedCases = failedCases + 1;
    end
    close all; diary off;
end

% --- Case 1B: static_sym_on ---
totalCases = totalCases + 1;
caseDir = fullfile(codeFolder, 'static_sym_on');
if exist(fullfile(caseDir, 'output.txt'), 'file')
    fprintf('[SKIP] TutorialShellsOptim / static_sym_on (already exists)\n');
    skippedCases = skippedCases + 1;
else
    fprintf('[RUN]  TutorialShellsOptim / static_sym_on ...\n');
    if ~exist(caseDir, 'dir'), mkdir(caseDir); end
    diary(fullfile(caseDir, 'output.txt'));
    try
        TutorialShellsOptim_(true);
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
%  2. TutorialShellsOptimDynamic — DYNAMIC, forceSymmetry = true
%% ========================================================================
codeFolder = fullfile(resultsRoot, 'TutorialShellsOptimDynamic');
if ~exist(codeFolder, 'dir'), mkdir(codeFolder); end

totalCases = totalCases + 1;
caseDir = fullfile(codeFolder, 'dynamic_sym_on');
if exist(fullfile(caseDir, 'output.txt'), 'file')
    fprintf('[SKIP] TutorialShellsOptimDynamic / dynamic_sym_on (already exists)\n');
    skippedCases = skippedCases + 1;
else
    fprintf('[RUN]  TutorialShellsOptimDynamic / dynamic_sym_on ...\n');
    if ~exist(caseDir, 'dir'), mkdir(caseDir); end
    diary(fullfile(caseDir, 'output.txt'));
    try
        TutorialShellsOptimDynamic_(true);
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
%  3. TutorialShellsOptimDeltaTheta — STATIC & DYNAMIC
%     deltaTheta0 = [0, -30, 15, 70]
%% ========================================================================
codeFolder = fullfile(resultsRoot, 'TutorialShellsOptimDeltaTheta');
if ~exist(codeFolder, 'dir'), mkdir(codeFolder); end

dThetaValues  = [0, -30, 15, 70];
dThetaLabels  = {'0', 'm30', '15', '70'};
optimCases_DT = {'STATIC', 'DYNAMIC'};
caseLabels_DT = {'static', 'dynamic'};

for ic = 1:length(optimCases_DT)
    for id = 1:length(dThetaValues)
        totalCases = totalCases + 1;
        caseName = sprintf('%s_dTheta_%s', caseLabels_DT{ic}, dThetaLabels{id});
        caseDir  = fullfile(codeFolder, caseName);

        if exist(fullfile(caseDir, 'output.txt'), 'file')
            fprintf('[SKIP] TutorialShellsOptimDeltaTheta / %s (already exists)\n', caseName);
            skippedCases = skippedCases + 1;
        else
            fprintf('[RUN]  TutorialShellsOptimDeltaTheta / %s ...\n', caseName);
            if ~exist(caseDir, 'dir'), mkdir(caseDir); end
            diary(fullfile(caseDir, 'output.txt'));
            try
                TutorialShellsOptimDeltaTheta_(optimCases_DT{ic}, dThetaValues(id));
                saveFigures(caseDir);
                casesRun = casesRun + 1;
            catch ME
                fprintf('ERROR: %s\n', ME.message);
                disp(getReport(ME));
                failedCases = failedCases + 1;
            end
            close all; diary off;
        end
    end
end

%% ========================================================================
%  4. DifferentThetasTutorialShellsOptim — STATIC & DYNAMIC
%     nOuter = [3, 5],  deltaTheta0 = [0, -30, 15, 70]
%% ========================================================================
codeFolder = fullfile(resultsRoot, 'DifferentThetasTutorialShellsOptim');
if ~exist(codeFolder, 'dir'), mkdir(codeFolder); end

nOuterValues   = [3, 5];
dThetaValues2  = [0, -30, 15, 70];
dThetaLabels2  = {'0', 'm30', '15', '70'};
optimCases_DT2 = {'STATIC', 'DYNAMIC'};
caseLabels_DT2 = {'static', 'dynamic'};

for io = 1:length(nOuterValues)
    for ic = 1:length(optimCases_DT2)
        for id = 1:length(dThetaValues2)
            totalCases = totalCases + 1;
            caseName = sprintf('%s_nOuter%d_dTheta_%s', ...
                caseLabels_DT2{ic}, nOuterValues(io), dThetaLabels2{id});
            caseDir  = fullfile(codeFolder, caseName);

            if exist(fullfile(caseDir, 'output.txt'), 'file')
                fprintf('[SKIP] DifferentThetasTutorialShellsOptim / %s (already exists)\n', caseName);
                skippedCases = skippedCases + 1;
            else
                fprintf('[RUN]  DifferentThetasTutorialShellsOptim / %s ...\n', caseName);
                if ~exist(caseDir, 'dir'), mkdir(caseDir); end
                diary(fullfile(caseDir, 'output.txt'));
                try
                    DifferentThetasTutorialShellsOptim_( ...
                        optimCases_DT2{ic}, nOuterValues(io), dThetaValues2(id));
                    saveFigures(caseDir);
                    casesRun = casesRun + 1;
                catch ME
                    fprintf('ERROR: %s\n', ME.message);
                    disp(getReport(ME));
                    failedCases = failedCases + 1;
                end
                close all; diary off;
            end
        end
    end
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
