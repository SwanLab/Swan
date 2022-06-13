
clc; clear variables; close all;
addpath(genpath(fileparts(mfilename('fullpath'))));

%% Main
fileName = 'test_cantileverEig'; % 'test_cantilever2';%'test_cantileverEig'; % 'test_cantilever'
settingsTopOpt = SettingsTopOptProblem(fileName);            

m = settingsTopOpt.problemData.femData.mesh;
x = m.coord(:,1);
y = m.coord(:,2);
xmax = max(x);

tol = 0.01;
a = 0;
%isNodeFixedRight = (abs(x -xmax)-0) < tol &  (abs(y - 0.5) - 0.1) < tol;
isNodeFixedDown  = (abs(x -a)-a) < tol &  (abs(y - 0) - 0.001) < tol;
isNodeFixedUp    = (abs(x -a)-a) < tol &  (abs(y - 1) - 0.001) < tol;
isNodeFixed = isNodeFixedUp | isNodeFixedDown;

settingsTopOpt.designVarSettings.isFixed.values = 1*ones(sum(isNodeFixed),1);
settingsTopOpt.designVarSettings.isFixed.nodes = find(isNodeFixed);

topOptProblem = TopOpt_Problem(settingsTopOpt);
topOptProblem.computeVariables;
topOptProblem.postProcess;
