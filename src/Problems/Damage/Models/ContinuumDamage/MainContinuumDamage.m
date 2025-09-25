clc,clear,close all

%% GENERAL SETTINGS
s.monitoring.set = false;
s.monitoring.print = false;

s.benchmark.mesh.type   = 'Rectangle';
s.benchmark.mesh.length = 1;
s.benchmark.mesh.width  = 1;
s.benchmark.mesh.lN     = 1;
s.benchmark.mesh.wN     = 1;

s.benchmark.bc.type   = 'DisplacementTractionY';
s.benchmark.bc.values = [0:1e-1:2];

s.matInfo.young   = 210;
s.matInfo.poisson = 0.3;

s.damageInfo.type = 'AT2';
s.damageInfo.r0 = 1e-10;
s.damageInfo.params.r1 = 1e10;
% s.damageInfo.params.hardening = -0.5;
% s.damageInfo.params.A = 1;
% s.damageInfo.params.qInf = 0;
s.damageInfo.params.w = 500;


s.tolerance = 1e-8;
s.maxIter   = 20;


%% RUN
tester = TestingContinuumDamage(s);
outputData = tester.compute();
outputData.inputParameters = s;

%% SAVE + PLOT
%save("~/Documents/GitHub/Swan/SENshear_SquareArea.mat","outputData") %ACTIVATE TO SAVE DATA!
ContinuumDamagePlotter(outputData);
