clc,clear,close all

s.monitoring.set = true;
s.monitoring.type = 'full'; %'reduced'
s.benchmark.type.mesh = '1Elem';
s.benchmark.N = 10;
s.benchmark.type.bc = 'displacementTraction';
s.benchmark.bcValues = [0:1e-4:1e-1];
s.matInfo.E  = 210;
s.matInfo.nu = 0.3;
s.matInfo.Gc = 5e-3;
s.matInfo.matType = 'PhaseFieldHomog';  %'PhaseFieldAnalytic'
s.matInfo.fileName = 'CircleMicroDamagePerimeter'; %'IsoMicroDamage','Circle/Square+MicroDamage+Area/Perimeter'
s.matInfo.degradation = 'PhaseFieldDegradation';
s.dissipInfo.type = 'PhaseFieldDissipationAT';
s.dissipInfo.pExp = 2;
s.l0 = 0.1;

tester = TestingPhaseField(s);
outputData = tester.compute();

PhaseFieldPlotter(outputData);

%save("X.mat","outputData")