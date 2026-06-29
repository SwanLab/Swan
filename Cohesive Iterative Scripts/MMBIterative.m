clear
clc
close all

%% Material

inputData.young   = 120e3;
inputData.poisson = 0.3;

%% Cohesive law

inputData.lawType = 'TractionBiliniarCoupled';

inputData.tau0Normal       = 15;
% inputData.tau0Shear        = 
inputData.firstCritEnergy  = 0.260;
inputData.secondCritEnergy = 1.002;
inputData.eta              = 2;

inputData.Kcoh = 1e6;

%% Boundary conditions
inputData.bcValues    = [0       :  5e-2    :  7,...
                             7       :   5e-2   :  8.2,...
                             8.2:         5e-2 : 10];


%% Problem definition

inputData.problemType = 'MMB';
inputData.solverType  = 'Newton';

%% Geometry

inputData.l = 150;
inputData.h = 1.55*2;

inputData.c = 63.18;

inputData.yCohLine    = 0.5*inputData.h;
inputData.xCohLineMax = 150-35;


%% Mesh
        a.fileName = 'CThickLever';
        s = FemDataContainer(a);
        baseMesh = s.mesh;
        inputData.baseMesh = baseMesh;

%% Cases to analyze

tau0ShearValues = [30 60];

results = struct();

%% Run simulations

for iCase = 1:length(tau0ShearValues)

    inputData.tau0Shear = tau0ShearValues(iCase);

    fprintf('\n=====================================\n');
    fprintf('Running case %d/%d: tau0Shear = %.1f MPa\n', ...
        iCase,length(tau0ShearValues),inputData.tau0Shear);
    fprintf('=====================================\n');

    problem = IterativeTutorialCohesive(inputData);

    results(iCase).tau0Shear = inputData.tau0Shear;
    results(iCase).problem   = problem;
    results(iCase).output    = problem.output;

end