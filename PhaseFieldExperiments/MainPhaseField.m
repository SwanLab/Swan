clc,clear,close all

s.monitoring.set = false;
s.monitoring.type = 'reduced'; %'reduced'
s.monitoring.print = false;
s.benchmark.type.mesh = 'SEN';%'1Elem';
s.benchmark.N = 10;
s.benchmark.type.bc = 'displacementShear';%'displacementTraction';
s.benchmark.bcValues = [0:1e-3:5e-3,5e-3:1e-3:2e-2];%[5e-3:0.1*1e-3:1e-2];%;%[0:1e-4:1e-1];
s.matInfo.E  = 210;
s.matInfo.nu = 0.3;
s.matInfo.Gc = 2.7e-3;
s.matInfo.matType = 'PhaseFieldAnalytic';  %'PhaseFieldAnalytic'
s.matInfo.fileName = 'IsoMicroDamage'; %'IsoMicroDamage','Circle/Square+MicroDamage+Area/Perimeter'
s.matInfo.degradation = 'PhaseFieldDegradation';
s.dissipInfo.type = 'PhaseFieldDissipationAT';
s.dissipInfo.pExp = 2;
s.l0 = 0.0075;

tester = TestingPhaseField(s);
outputData = tester.compute();

PhaseFieldPlotter(outputData);

save("SENshear1e-3.mat","outputData")