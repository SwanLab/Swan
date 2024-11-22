clc,clear,close all

% load("SEN1e-3_0.0015_1e-12.mat")
% s.initialGuess.u = outputData.displacement.field;
% s.initialGuess.phi = outputData.damage.field;

s.monitoring.set = false;
s.monitoring.type = 'full'; %'reduced'
s.monitoring.print = true;
s.benchmark.type.mesh = 'SEN';%'1Elem';
s.benchmark.N = 10;
s.benchmark.type.bc = 'displacementShear';%'displacementTraction';

%% TEST CASES 
% SEN CASE 1
s.benchmark.bcValues = [0:1e-4:8e-3,8e-3:1e-5:1e-1];
s.matInfo.E  = 210;
s.matInfo.nu = 0.3;
s.matInfo.Gc = 2.7e-3;
s.l0 = 0.0015;

% % SEN CASE 2
% s.benchmark.bcValues = [0:1e-3:8e-3,8e-3:1e-3:4e-1];
% s.matInfo.E  = 210;
% s.matInfo.nu = 0.3;
% s.matInfo.Gc = 2.7e-3;
% s.l0 = 0.0015;
% CANVIAR SAVE PER: save("SEN1e-3_0.0015_1e-12.mat","outputData")

% % SEN CASE 3
% s.benchmark.bcValues = [0:1e-4:8e-3,8e-3:1e-5:3e-1];
% s.matInfo.E  = 210;
% s.matInfo.nu = 0.3;
% s.matInfo.Gc = 2.7e-3;
% s.l0 = 0.01;
% CANVIAR SAVE PER: save("SEN1e-5_0.01_1e-12.mat","outputData")

s.matInfo.matType = 'PhaseFieldAnalytic';%'PhaseFieldHomog';  %'PhaseFieldAnalytic'
s.matInfo.fileName = 'CircleMicroDamagePerimeter'; %'IsoMicroDamage','Circle/Square+MicroDamage+Area/Perimeter'
s.matInfo.degradation = 'PhaseFieldDegradation';
s.dissipInfo.type = 'PhaseFieldDissipationAT';
s.dissipInfo.pExp = 2;


tester = TestingPhaseField(s);
outputData = tester.compute();
save("SEN1e-5_0.0015_1e-12.mat","outputData") %ACTIVATE TO SAVE DATA!

PhaseFieldPlotter(outputData);



%% NO TOCAR
% Lshape
% s.benchmark.bcValues = [0.01:0.01:0.35,0.32:-0.03:0.02];%[5e-3:0.1*1e-3:1e-2];%;%[0:1e-4:1e-1];
% s.matInfo.E  = 25.85*10e3;
% s.matInfo.nu = 0.18;
% s.matInfo.Gc = 95*10e-4;
% s.l0 = 5;
