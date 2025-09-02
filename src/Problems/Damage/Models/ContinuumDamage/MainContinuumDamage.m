clc,clear,close all

%% GENERAL SETTINGS
s.monitoring.set = true;
s.monitoring.print = true;

s.benchmark.mesh.length = 100;
s.benchmark.mesh.width  = 1;
s.benchmark.mesh.lN     = 100;
s.benchmark.mesh.wN     = 100;
s.benchmark.mesh.type   = 'Sample';

s.benchmark.bc.type     = 'DisplacementTractionX';
s.benchmark.bc.bcValues = [0:1e-2:1e-1,1e-1:1e-3:2e-1];

s.matInfo.young   = 210;
s.matInfo.poisson = 0.3;

s.damageInfo.type = 'Linear';
s.damageInfo.r0 = 1;
s.damageInfo.params.r1 = 3;
s.damageInfo.params.hardening = -0.5;
% obj.damageInfo.params.A = 0.1;
% obj.damageInfo.params.w = 0.1 ;
%obj.damageInfo.params.qInf = 30;

s.tolerance = 1e-6;
s.maxIter   = 20;


%% RUN
tester = TestingContinuumDamage(s);
outputData = tester.compute();
outputData.inputParameters = s;

%% SAVE + PLOT
%save("~/Documents/GitHub/Swan/SENshear_SquareArea.mat","outputData") %ACTIVATE TO SAVE DATA!
%ContinuumDamagePlotter(outputData);
