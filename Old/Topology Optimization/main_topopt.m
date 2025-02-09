
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
%fileName = 'CantileverTriangle_Case_1_1_1';
%fileName = 'CantileverTetrahedra_Case_1_1_1';
%fileName = 'InteriorPerimeter';
fileName = 'MinPoissLevelSet';
% settings = Settings('Case_RefactoringSettings_OLD');
% settingsTopOpt = SettingsTopOptProblem('Case_RefactoringSettings_A',settings);
% settings = Settings('Case_RefactoringSettingsMICRO_OLD');
% settingsTopOpt = SettingsTopOptProblem('Case_RefactoringSettingsMICRO_A',settings);

% settings = Settings('Case_RefactoringSettingsMICRO_OLD_2');
% settingsTopOpt = SettingsTopOptProblem('CaseBenchmark_JSON_A.json',settings);
% 
% settings = Settings('Case_RefactoringSettings_OLD');
% settingsTopOpt = SettingsTopOptProblem('CaseBenchmark_JSON_B.json',settings);
%settings = Settings('CantileverTriangle_Case_1_1_1InteriorPerimeter');
%settings = Settings('InteriorPerimeter');
%settings = Settings('CantileverHexahedra_Case_1_1_1');
%settings = Settings('MinPoissLevelSet');
 %settings = Settings('CantileverTetrahedraFine_Case_1_1_1');
%settings = Settings('CantileverHexahedraCoarse_Case_1_1_1');
%settings = Settings('ImprovedBridgeSYM_Case_1_1_1');
%settings.printing = false;
%settings.plotting = false;
%settings.monitoring = false;
settingsTopOpt = SettingsTopOptProblem(fileName);            

% settingsTopOpt = SettingsTopOptProblem('CaseBenchmark_JSON_B.json');


topOptProblem = TopOpt_Problem(settingsTopOpt);
topOptProblem.computeVariables;
topOptProblem.postProcess;

%close all
