
clc; clear variables; close all;
addpath(genpath(fileparts(mfilename('fullpath'))));

%% Main
fileName = 'test_cantileverEig'; % 'test_cantilever2';%'test_cantileverEig'; % 'test_cantilever'
settingsTopOpt = SettingsTopOptProblem(fileName);            
topOptProblem = TopOpt_Problem(settingsTopOpt);
topOptProblem.computeVariables;
topOptProblem.postProcess;
