%% Swan - Wiki
% The student's guide to clean code development
% Task 8: UML of Topology Optimization code in Swan repository

% Instructions: run the following code, selecting previously the 'Swan'
% main folder as your current matlab path

% fclose('all');
% rmdir('Output','s')
clc
close all

%%

tic

s.testName = 'test_micro3d';
t = TopOptComputer(s);
t.compute();

toc