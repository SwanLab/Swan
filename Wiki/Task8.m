%% Swan - Wiki
% The student's guide to clean code development
% Task 8: UML of Topology Optimization code in Swan repository

% Instructions: run the following code, selecting previously the 'Swan'
% main folder as your current matlab path

close all
clear

s.testName = 'test_cantilever2';
t = TopOptComputer(s);
t.compute();

s.filename = s.testName;
s.type     = 'GiD';
s.mesh     = t.computation.designVariable.mesh;
s.fValues  = t.computation.designVariable.value;
Result     = P1Function(s);
Result.print(s);