% clear
% clc
% close all

%% Material

inputData.young   = 120e9;
inputData.poisson = 0.3;

% E  = inputData.young;
% nu = inputData.poisson;
% 
% inputData.young   = E/(1-nu^2);
% inputData.poisson = nu/(1-nu);


%% Cohesive law

inputData.lawType = 'TractionBiliniarCoupled';

% Parameters that remain constant
inputData.tau0Normal        = 15e6;
inputData.firstCritEnergy   = 260;
inputData.secondCritEnergy  = 1002;
inputData.eta               = 2;

% Alternative parameters
inputData.jumpCrit  = 1.25e-7;
inputData.jumpFinal = 0.025e-3;

inputData.Kcoh = 1e13;

%% Boundary conditions


% % tau3o 30
bcValues{1} = [ ...
    0     : 5e-4     :6.3e-3,...
    6.3e-3  : 1e-5    :7.3e-3,...
    7e-3  : 5e-5     :9e-3];

% % tau3o 60
bcValues{2} = [ ...
    0     : 5e-4     :6.8e-3,...
    6.8e-3  : 1e-5    :7.7e-3,...
    7.7e-3  : 5e-5     :9e-3];

% tau3o 120
bcValues{3} = [ ...
    0     : 5e-4     :6.2e-3,...
    6.2e-3: 1e-4 : 7.75e-3,...
    7.75e-3  : 1e-5    :8.3e-3,...
    8.3e-3  :   5e-5   :9e-3];



%% Problem definition

inputData.problemType = 'EndNotchedFlex';
inputData.solverType  = 'Newton';

%% Geometry

inputData.l = 150e-3;
inputData.h = 1.55*2e-3;

inputData.yCohLine    = 0.5*inputData.h;
inputData.xCohLineMax = 115e-3;

%% Mesh
    inputData.nx = 1500;
    inputData.ny = 10;

%% Cases to analyze


tau0ShearValues = [30e6, 60e6, 120e6];

results = struct();

%% Run simulations

for iCase = 1:length(tau0ShearValues)
    inputData.bcValues = bcValues{iCase};
    inputData.tau0Shear = tau0ShearValues(iCase);

    fprintf('\n=====================================\n');
    fprintf('Running case %d/%d: tau0Shear = %.1f MPa\n', ...
        iCase, length(tau0ShearValues), inputData.tau0Shear/1e6);
    fprintf('=====================================\n');

    problem = IterativeTutorialCohesive(inputData);

    results(iCase).tau0Shear = inputData.tau0Shear;
    results(iCase).problem    = problem;
    results(iCase).output     = problem.output;

end

