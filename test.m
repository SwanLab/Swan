clear all; close all; clc
%% TEST
% - 

% Load the results for 2-d test
load './test-2d/d_u.mat'

% Parent directory
[parentdir,~,~] = fileparts(pwd);

% Run Main.m
obj = Main;

if sum(abs(obj.variables.displacement - d_u)) < 1e-6
    disp('TEST PASSED');
else
    disp('TEST FAILED');
end