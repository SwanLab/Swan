clc,clear,close all

s.monitoring.set = true;
s.monitoring.type = 'full'; %'reduced'
s.monitoring.print = true;
s.benchmark.type.mesh = 'Lshape';%'1Elem';
s.benchmark.N = 10;
s.benchmark.type.bc = 'Lshape';%'displacementTraction';
s.benchmark.bcValues = [0.01:0.01:0.35,0.32:-0.03:0.02];%[0:1e-4:5e-3,5e-3:1e-4:2e-1];%[5e-3:0.1*1e-3:1e-2];%;%[0:1e-4:1e-1];
s.matInfo.E  = 25.85*10e3;
s.matInfo.nu = 0.18;
s.matInfo.Gc = 95*10e-6;
s.matInfo.matType = 'PhaseFieldAnalytic';%'PhaseFieldHomog';  %'PhaseFieldAnalytic'
s.matInfo.fileName = 'CircleMicroDamagePerimeter'; %'IsoMicroDamage','Circle/Square+MicroDamage+Area/Perimeter'
s.matInfo.degradation = 'PhaseFieldDegradation';
s.dissipInfo.type = 'PhaseFieldDissipationAT';
s.dissipInfo.pExp = 2;
s.l0 = 5;

tester = TestingPhaseField(s);
outputData = tester.compute();
save("SEN1e-5_0.0015_1e-12.mat","outputData") %ACTIVATE TO SAVE DATA!

PhaseFieldPlotter(outputData);



