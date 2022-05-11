
clc; clear variables; close all;
addpath(genpath(fileparts(mfilename('fullpath'))));

%% Main
fileName = 'test_cantilever2';
settingsTopOpt = SettingsTopOptProblem(fileName);            

topOptProblem = TopOpt_Problem(settingsTopOpt);
topOptProblem.computeVariables;
topOptProblem.postProcess;