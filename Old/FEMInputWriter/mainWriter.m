%% Main Writer - Cantilever Test (2)

clear;
clc;

% Data input
s.testName = 'jaCantilever.m';
s.x1       = 2;
s.y1       = 1;
s.N        = 50;
s.M        = 25;
s.P        = -100;
s.DoF      = 2;

% File writing
FEMWriter = FEMInputWriter(s);
FEMWriter.createTest;