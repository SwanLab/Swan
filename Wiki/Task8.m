%% Swan - Wiki
% The student's guide to clean code development
% Task 8: UML of Topology Optimization code in Swan repository

% Instructions: run the following code, selecting previously the 'Swan'
% main folder as your current matlab path
clc
clear all
close all
s.testName = 'test_grip';
s.filename = s.testName;
            scaleParameters.eta = 0.25;
            scaleParameters.beta = 1.5;
            scaleParameters.rad = 0.28;
            save('scaleParameters.mat', 'scaleParameters')
t = TopOptComputer(s);
t.compute();

p.mesh    = t.computation.designVariable.mesh;
p.fValues = t.computation.designVariable.value;
Result = P1Function(p);
Result.print(s);