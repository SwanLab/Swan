clc,clear,close all

%% INITIAL CONDITIONS FOR CONTINUING AN UNFINISHED ANALYSIS
% load("SENshearAT2_v2.mat")
% s.initialGuess.u = outputData.displacement.field;
% s.initialGuess.phi = outputData.damage.field.fun;

%% GENERAL SETTINGS

s.monitoring.set = true;
s.monitoring.type = 'full'; %'reduced'
s.monitoring.print = true;

s.tolerance.u = 1e-6;
s.tolerance.phi = 1e-6;
s.tolerance.stag = 1e-6;
s.maxIter.u = 100;
s.maxIter.phi = 300;
s.maxIter.stag = 300;

s.benchmark.mesh.type   = '1Elem';
% s.benchmark.mesh.length = 200;
% s.benchmark.mesh.width  = 10;
% s.benchmark.mesh.lN     = 200;
% s.benchmark.mesh.wN     = 10;

s.benchmark.bc.phi.type = 'DamageFree';
s.benchmark.bc.phi.isAdaptive = false;

s.benchmark.bc.u.type   = 'DisplacementTractionYClamped';
s.benchmark.bc.u.isAdaptive = false;
s.benchmark.bc.u.values =  [0:1e-3:0.1];
% s.benchmark.bc.u.isAdaptive = true;
% s.benchmark.bc.u.initialValue = 0;
% s.benchmark.bc.u.finalValue = 0.03;

s.matInfo.matType = 'Analytic'; %'Analytic','Homogenized'
s.matInfo.degradationType = 'PhaseField'; %'PhaseField','SIMPALL'
s.matInfo.degradationSubType = 'General'; %'AT','AT2linear','General'
s.matInfo.fileName = 'HorizontalCrack'; 
s.matInfo.young   = 210;
s.matInfo.poisson = 0.3;
s.matInfo.Gc      = 5e-3;
s.matInfo.sigmaMax = 1.984;
s.l0 = 0.1;
s.matInfo.params.coeffs = [(4/pi)*(s.matInfo.Gc*s.matInfo.young)/(s.matInfo.sigmaMax^2 * s.l0), -0.5]; %(4/pi)
s.matInfo.params.exp = 2;

s.dissipInfo.type = 'FullQuadratic';
s.dissipInfo.constant = 8/3; % 8/5 AT3 / 2 AT2 / 8/3 AT1 / pi Rational 
s.dissipInfo.pExp = 1;
s.dissipInfo.xi = 2; % 0 AT2 / 1 AT1 / 2 Mix in type FullQuadratic
s.solver.type = 'Gradient';
s.solver.tau  = 150;


%% Orientation plots
angleInfo.alpha = [];
angleInfo.beta  = [];
angleInfo.theta = [];
angleInfo.chi   = [];
angleInfo.phi   = [];
angleInfo.phiOld= [];
save('angleInfo','angleInfo')

%% RUN
tester = TestingPhaseField(s);
outputData = tester.compute();
outputData.inputParameters = s;

%% SAVE + PLOT
%save("1ElemRational1_v2.mat",'outputData') %ACTIVATE TO SAVE DATA!
PhaseFieldPlotter(outputData);

%% ORIENTATION PLOTS

for i=1:4

angleInfo = outputData.angleInfo;
aData = angleInfo.alpha(i,:);
bData = angleInfo.beta(i,:);
cData = angleInfo.chi(i,:);
tData = angleInfo.theta(i,:);

% bData = unwrap(2*bData)/2;
tData = wrapTo2Pi(unwrap(2*tData)/2);
%aData = wrapTo2Pi(unwrap(2*(bData+tData))/2);
aData = bData+tData;
%bData(abs(bData - 180) < 1e-8) = 0;
%tData(abs(tData - 180) < 1e-8) = 0;

figure(100+i)
hold on
yyaxis right
p1 = plot(angleInfo.phi(i,:),'Color','Y');
p2 = plot(angleInfo.phiOld(i,:));
ylim([-0.1 1.1])
legend('phi')
yyaxis left
p3 = plot(rad2deg(aData),'Color','R');
p4 = plot(rad2deg(bData),'Color','B');
p5 = plot(rad2deg(tData),'Color','G');
p6 = plot(rad2deg(cData),'Color','K');
ylim([-5 365])
legend([p1,p2,p3,p4,p5,p6],'phi','phiOld','alpha','beta','theta','chi')
end



