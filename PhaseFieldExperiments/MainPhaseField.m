clc,clear,close all

s.monitoring.set = false;
s.monitoring.type = 'full'; %'reduced'
s.monitoring.print = true;
s.benchmark.type.mesh = 'SEN';%'1Elem';
s.benchmark.N = 10;
s.benchmark.type.bc = 'displacementShear';%'displacementTraction';
s.benchmark.bcValues = [0:1e-4:3e-2,3e-2:1e-5:3.5e-2];%[0:1e-4:5e-3,5e-3:1e-4:2e-1];%[5e-3:0.1*1e-3:1e-2];%;%[0:1e-4:1e-1];
s.matInfo.E  = 210;
s.matInfo.nu = 0.3;
s.matInfo.Gc = 2.7e-3;
s.matInfo.matType = 'PhaseFieldAnalytic';%'PhaseFieldHomog';  %'PhaseFieldAnalytic'
s.matInfo.fileName = 'CircleMicroDamagePerimeter'; %'IsoMicroDamage','Circle/Square+MicroDamage+Area/Perimeter'
s.matInfo.degradation = 'PhaseFieldDegradation';
s.dissipInfo.type = 'PhaseFieldDissipationAT';
s.dissipInfo.pExp = 2;
s.l0 = 0.0015;

tester = TestingPhaseField(s);
outputData = tester.compute();

PhaseFieldPlotter(outputData);

%save("SEN1e-4_0.0075_1e-8.mat","outputData") %ACTIVATE TO SAVE DATA!

