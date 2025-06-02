% Stoke Problem Optimization
% clc
% clear;
close all;

% Read Network data
%load("StokesNetwork.mat");
load("StokesNetworkE0.5e5N36HL6MaxAoA12.88.mat");

% Initialization
m   = 0.09;
p   = 0.8;
t   = 0.4;
AoA = 1;

d.features      = [m,p,t,AoA];
d.learningRate  = 0.015; %0.06  %1.5 mínim per LE fixe
d.optimizer     = opt;
d.tol           = 1e-6;
d.lowerBC = [0.0, 0.0, 0.1, -2.5]; % 0.0, 0.0, 0.1, -2.5
d.upperBC = [0.09, 0.9, 0.4, 15.0];

% Compute the Optimization
stokeOpt = AirfoilOptimizer(d);
stokeOpt.computeOptAirfoilParams();close 

% Results
OptimalParams = stokeOpt.optimalParams;

% Plot the Optimal Airfoil
% figure;
% AirfoilOptimizer.plotAirfoilContour(OptimalParams,0)
%stokeOpt.plotEEvolution();
%% Plot E vs iterations

%stokeOpt.plotEEvolution();

%% Generate the Airfoil Shape Optimization Video

stokeOpt.generateAFSOPVideo();
 close all
stokeOpt.generateVelVideo();
 close all
stokeOpt.generatePVideo();
 close all

system('shutdown /s /t 60');
