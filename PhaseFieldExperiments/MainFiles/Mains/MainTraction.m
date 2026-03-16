clc,clear,close all

%% GENERAL SETTINGS

s.monitoring.set = false;
s.monitoring.type = 'full'; %'reduced'
s.monitoring.print = false;

s.tolerance.u = 1e-6;
s.tolerance.phi = 1e-6;
s.tolerance.stag = 1e-6;
s.maxIter.u = 100;
s.maxIter.phi = 300;
s.maxIter.stag = 300;

s.benchmark.mesh.type   = 'SENtraction';
s.benchmark.bc.u.type   = 'DisplacementTractionYClamped';
s.benchmark.bc.u.values =  [0:1e-5:0.005,0.005:1e-6:0.0065];
s.benchmark.bc.phi.type = 'DamageFree';

s.matInfo.matType = 'Analytic'; %'Analytic','Homogenized'
s.matInfo.degradationType = 'PhaseField'; %'PhaseField','SIMPALL'
s.matInfo.degradationSubType = 'AT'; %'AT','ATSplit','Rational','General'
s.matInfo.fileName = 'HoneycombBenchmark02'; 
s.matInfo.young   = 210;
s.matInfo.poisson = 0.3;
s.matInfo.Gc      = 2.7e-3;
s.matInfo.sigmaMax = 2.44542;
s.l0 = 0.01;
s.matInfo.params.coeffs = [(4/pi)*(s.matInfo.Gc*s.matInfo.young)/(s.matInfo.sigmaMax^2 * s.l0), -0.5]; %(4/pi)
s.matInfo.params.exp = 2;

s.dissipInfo.type = 'AT';
s.dissipInfo.constant = 8/3; % 2 AT2 / 8/3 AT1 / pi Rational 
s.dissipInfo.pExp = 1;
s.dissipInfo.xi = 2; % 0 AT2 / 1 AT1 / 2 Mix in type FullQuadratic
s.solver.type = 'Gradient';
s.solver.tau  = 150;


%% RUN
% tester = TestingPhaseField(s);
% outputData = tester.compute();
% outputData.inputParameters = s;
% save("SENtractionNewMeshAT1.mat",'outputData')
% 
% s.matInfo.degradationSubType = 'General';
% s.matInfo.sigmaMax = 2.44542;
% s.matInfo.params.coeffs = [(4/pi)*(s.matInfo.Gc*s.matInfo.young)/(s.matInfo.sigmaMax^2 * s.l0), -0.5]; 
% s.dissipInfo.constant = pi;
% tester = TestingPhaseField(s);
% outputData = tester.compute();
% outputData.inputParameters = s;
% save("SENtractionNewMeshRational24.mat",'outputData')

s.matInfo.degradationSubType = 'General'; 
s.matInfo.sigmaMax = 5;
s.matInfo.params.coeffs = [(4/pi)*(s.matInfo.Gc*s.matInfo.young)/(s.matInfo.sigmaMax^2 * s.l0), -0.5]; 
s.dissipInfo.constant = pi;
tester = TestingPhaseField(s);
outputData = tester.compute();
outputData.inputParameters = s;
save("SENtractionNewMeshRational5.mat",'outputData')

% s.matInfo.matType = 'Homogenized'; %'Analytic','Homogenized'
% s.matInfo.fileName = 'HexagonBenchmark03';
% s.dissipInfo.constant = 8/3;
% tester = TestingPhaseField(s);
% outputData = tester.compute();
% outputData.inputParameters = s;
% save("SENtractionNewMeshHexagon.mat",'outputData')
% 
% s.matInfo.matType = 'Homogenized'; %'Analytic','Homogenized'
% s.matInfo.fileName = 'HoneycombBenchmark03';
% s.dissipInfo.constant = 8/3;
% tester = TestingPhaseField(s);
% outputData = tester.compute();
% outputData.inputParameters = s;
% save("SENtractionNewMeshHoneycomb.mat",'outputData')