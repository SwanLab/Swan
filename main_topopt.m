
clc; clear variables; close all;
addpath(genpath(fileparts(mfilename('fullpath'))));

%% Test
% run('PlottingTests.m')
% run('FemTests.m');
% run('TopOptTests.m');
% run('UnfittedIntegrationTests.m')
% run('ExploringSettingsTests.m')
% run('AllTests.m')
clear variables;

%% Main
% settings = Settings('Case_RefactoringSettings_OLD');
% settingsTopOpt = SettingsTopOptProblem('Case_RefactoringSettings_A',settings);

% settings = Settings('Case_RefactoringSettingsMICRO_OLD');
% settingsTopOpt = SettingsTopOptProblem('Case_RefactoringSettingsMICRO_A',settings);

% settings = Settings('Case_RefactoringSettingsMICRO_OLD_2');
% settingsTopOpt = SettingsTopOptProblem('CaseBenchmark_JSON_A.json',settings);

settings = Settings('Case_RefactoringSettings_OLD');
settingsTopOpt = SettingsTopOptProblem('CaseBenchmark_JSON_B.json',settings);

% settings = Settings('CantileverTetrahedra_Case_1_1_1');
% SettingsTranslator().translate(settings);
% settingsTopOpt = SettingsTopOptProblem('CantileverTetrahedra_Case_1_1_1.json');

topOptProblem = TopOpt_Problem(settingsTopOpt);
topOptProblem.computeVariables;
topOptProblem.postProcess;

close all
