clc,clear,close all

%% INITIAL CONDITIONS
% load("SEN1e-5_0.0015_1e-8v3.mat")
% s.initialGuess.u = outputData.displacement.field;
% s.initialGuess.phi = outputData.damage.field;

%% SETTINGS
s.monitoring.set = true;
s.monitoring.type = 'full'; %'reduced'
s.monitoring.print = true;
s.benchmark.N = 10;
s.tolerance.u = 1e-13;
s.tolerance.phi = 1e-8;
s.tolerance.stag = 1e-8;

% 1Elem
% s.benchmark.type.mesh = '1Elem';
% s.benchmark.type.bc = 'displacementTraction';
% s.benchmark.bcValues = [1e-10:1e-3:1e-1];
% s.matInfo.E  = 210;
% s.matInfo.nu = 0.3;
% s.matInfo.Gc = 5e-3;
% s.l0 = 0.1;

% SEN Traction
s.benchmark.type.mesh = 'SENtraction';
s.benchmark.type.bc = 'displacementTraction';
%s.benchmark.bcValues = [1e-4:1e-4:5e-3,5e-3:1e-5:6e-3]; %AT2
s.benchmark.bcValues = [1e-4:1e-3:1e-2,1e-2:1e-4:1e-1];
%s.benchmark.bcValues = [1e-3:1e-4:1.5e-2]; %AT1
s.matInfo.E  = 210;
s.matInfo.nu = 0.3;
s.matInfo.Gc = 2.7e-3;
s.l0 = 0.0015; %0.0015;

% % SEN Shear
% s.benchmark.type.mesh = 'SENshear';
% s.benchmark.type.bc = 'displacementShear';
% s.benchmark.bcValues = [0:1e-3:1e-1,1e-1:1e-3:1.2e-1];
% s.matInfo.E  = 210;
% s.matInfo.nu = 0.3;
% s.matInfo.Gc = 2.7e-3;
% s.l0 = 0.01;

% SEN Mixed
% s.benchmark.type.mesh = 'SENmixed';
% s.benchmark.type.bc = 'displacementMixed';
% %s.benchmark.bcValues = [1e-3:1e-3:1.5e-2]; %AT2
% s.benchmark.bcValues = [1e-3:1e-3:1.5e-2]; %AT1
% s.benchmark.bcValues = [1e-4:1e-3:1e-2,1e-2:1e-4:1e-1]; %HomogPeri
% s.matInfo.E  = 210;
% s.matInfo.nu = 0.3;
% s.matInfo.Gc = 2.7e-3;
% s.l0 = 0.0015; %0.01

% Lshape
% s.benchmark.type.mesh = 'Lshape';
% s.benchmark.type.bc = 'Lshape';
% s.benchmark.bcValues = [0.01:0.01:0.35,0.32:-0.03:0.02];
% s.matInfo.E  = 25.85*10e3;
% s.matInfo.nu = 0.18;
% s.matInfo.Gc = 95*10e-6*100;
% s.l0 = 5;

s.matInfo.matType = 'PhaseFieldHomog';%'PhaseFieldHomog';  %'PhaseFieldAnalytic'
s.matInfo.fileName = 'CircleMicroDamageArea'; %'IsoMicroDamage','Circle/Square+MicroDamage+Area/Perimeter'
s.matInfo.degradation = 'PhaseFieldDegradation';
s.dissipInfo.type = 'PhaseFieldDissipationAT';
s.dissipInfo.pExp = 2;
s.solverType = 'Gradient'; %'Newton'

%% RUN
tester = TestingPhaseField(s);
outputData = tester.compute();
outputData.inputParameters = s;
save("/home/gerard/Documents/GitHub/Swan/PhaseFieldExperiments/ResultsOctober/SENtraction/" + ...
      "SENtraction_CircleArea_Gradient2.mat","outputData") %ACTIVATE TO SAVE DATA!

PhaseFieldPlotter(outputData);


