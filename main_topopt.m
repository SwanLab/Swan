
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
settings = Settings('Case_RefactoringSettings_OLD');
settingsTopOpt = SettingsTopOptProblem('Case_RefactoringSettings_A',settings);

topOptProblem = TopOpt_Problem(settingsTopOpt);
topOptProblem.computeVariables;
topOptProblem.postProcess;

close all
