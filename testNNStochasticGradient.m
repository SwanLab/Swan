clear;
clc;

rng(0);

%% ============================================================
%  LOAD DATA
% ============================================================

if ~exist('HomogNN.mat','file')

    error( ...
        ['HomogNN.mat was not found. ', ...
         'The test only needs the stored data.']);

end

S = load('HomogNN.mat','data');

data = ...
    S.data;


%% ============================================================
%  CREATE A FRESH NETWORK
%
%  No training is performed.
%
%  We use TANH first because it is smooth and therefore gives
%  a cleaner finite-difference verification.
% ============================================================

networkParams.hiddenLayers = ...
    [150 200 300 200 150 50];

networkParams.HUtype = ...
    'tanh';

networkParams.OUtype = ...
    'linear';

networkParams.data = ...
    data;


network = ...
    Network(networkParams);


%% ============================================================
%  LEARNABLE VARIABLES
% ============================================================

thetaVar = ...
    network.getLearnableVariables();

theta0 = ...
    thetaVar.thetavec;


fprintf('\n');
fprintf('=============================================================\n');
fprintf(' NN STOCHASTIC GRADIENT CHECK\n');
fprintf('=============================================================\n');

fprintf( ...
    'Number of parameters: %d\n', ...
    numel(theta0));

fprintf( ...
    'Number of training samples: %d\n', ...
    size(data.Xtrain,1));

fprintf( ...
    'Number of network inputs: %d\n', ...
    size(data.Xtrain,2));

fprintf( ...
    'Number of outputs: %d\n', ...
    size(data.Ytrain,2));


%% ============================================================
%  LOSS FUNCTIONAL
% ============================================================

sLoss.costType = ...
    'L2';

sLoss.designVariable = ...
    thetaVar;

sLoss.network = ...
    network;

sLoss.data = ...
    data;


loss = ...
    LossFunctional(sLoss);


%% ============================================================
%  L2 REGULARIZATION
% ============================================================

sReg.designVariable = ...
    thetaVar;


regularization = ...
    Sh_Func_L2norm(sReg);


%
% Use a non-zero lambda so that the stochastic CostNN path
% also exercises the regularization contribution.
%
lambda = ...
    1e-4;


%% ============================================================
%  TOTAL COST
% ============================================================

sCost.shapeFunctions = ...
    {loss,regularization};

sCost.weights = ...
    [1,lambda];


cost = ...
    CostNN(sCost);


%% ============================================================
%  FREEZE THE MINI-BATCH
%
%  This is essential.
%
%  moveBatch = false means that every evaluation below must
%  use the same mini-batch.
% ============================================================

cost.setBatchMover(false);


%% ============================================================
%  FIRST STOCHASTIC EVALUATION
% ============================================================

cost.computeStochasticFunctionAndGradient(theta0);

J0 = ...
    cost.value;

g = ...
    cost.gradient;


if ~isequal(size(g),size(theta0))

    error( ...
        'Gradient shape mismatch: theta is %s, gradient is %s.', ...
        mat2str(size(theta0)), ...
        mat2str(size(g)));

end


fprintf('\nFixed mini-batch evaluation:\n');

fprintf( ...
    'J_B(theta) = %.16e\n', ...
    J0);

fprintf( ...
    '||grad J_B||_2 = %.16e\n', ...
    norm(g(:),2));


%% ============================================================
%  VERIFY THAT moveBatch = false REALLY KEEPS THE SAME BATCH
% ============================================================

cost.setBatchMover(false);

cost.computeStochasticFunctionAndGradient(theta0);

Jrepeat = ...
    cost.value;


batchRepeatError = ...
    abs(Jrepeat-J0);


fprintf('\nRepeated evaluation at the same theta:\n');

fprintf( ...
    '|J_B(theta) - J_B(theta) repeated| = %.16e\n', ...
    batchRepeatError);


if batchRepeatError > 1e-14

    warning( ...
        ['Repeated stochastic evaluation changed even though ', ...
         'moveBatch=false. The batch may not be fixed.']);

end


%% ============================================================
%  RANDOM DIRECTION
% ============================================================

v = ...
    randn(size(theta0));

v = ...
    v / norm(v(:),2);


%
% Analytical directional derivative:
%
%       D_v J_B = grad(J_B)^T v
%
dJanalytical = ...
    sum(g(:).*v(:));


fprintf('\nAnalytical directional derivative:\n');

fprintf( ...
    'grad(J_B)^T v = %.16e\n', ...
    dJanalytical);


%% ============================================================
%  CENTRAL FINITE DIFFERENCE
%
%                  J_B(theta+h*v) - J_B(theta-h*v)
%       D_v J_B ~= -----------------------------------
%                                2h
%
%  All evaluations must use THE SAME mini-batch.
% ============================================================

hValues = ...
    [1e-2 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8];


absErrors = ...
    zeros(size(hValues));

relErrors = ...
    zeros(size(hValues));


fprintf('\n');
fprintf('--------------------------------------------------------------------------\n');
fprintf('       h              FD derivative          abs error          rel error\n');
fprintf('--------------------------------------------------------------------------\n');


for ih = 1:numel(hValues)

    h = ...
        hValues(ih);


    %% --------------------------------------------------------
    % theta + h*v
    % ---------------------------------------------------------

    thetaPlus = ...
        theta0 + h*v;

    cost.setBatchMover(false);

    cost.computeStochasticFunctionAndGradient( ...
        thetaPlus);

    Jplus = ...
        cost.value;


    %% --------------------------------------------------------
    % theta - h*v
    % ---------------------------------------------------------

    thetaMinus = ...
        theta0 - h*v;

    cost.setBatchMover(false);

    cost.computeStochasticFunctionAndGradient( ...
        thetaMinus);

    Jminus = ...
        cost.value;


    %% --------------------------------------------------------
    % Central finite difference
    % ---------------------------------------------------------

    dJFD = ...
        (Jplus-Jminus)/(2*h);


    %% --------------------------------------------------------
    % Errors
    % ---------------------------------------------------------

    absErrors(ih) = ...
        abs(dJFD-dJanalytical);


    relErrors(ih) = ...
        absErrors(ih) ...
        / ...
        max(abs(dJanalytical),eps);


    fprintf( ...
        '%12.1e    %20.12e    %12.4e    %12.4e\n', ...
        h, ...
        dJFD, ...
        absErrors(ih), ...
        relErrors(ih));

end


%% ============================================================
%  RESTORE ORIGINAL PARAMETERS
% ============================================================

thetaVar.thetavec = ...
    theta0;


%% ============================================================
%  RESULT
% ============================================================

[minRelError,idBest] = ...
    min(relErrors);


fprintf('--------------------------------------------------------------------------\n');

fprintf( ...
    '\nBest relative error = %.6e at h = %.1e\n', ...
    minRelError, ...
    hValues(idBest));


if minRelError < 1e-5

    fprintf('\nRESULT: PASS\n');

elseif minRelError < 1e-3

    fprintf( ...
        '\nRESULT: ACCEPTABLE, but inspect convergence with h.\n');

else

    fprintf( ...
        '\nRESULT: FAIL - stochastic gradient must be investigated.\n');

end


fprintf('=============================================================\n');