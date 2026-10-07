function testSGD()

clc;

nIn  = 2;
nOut = 6;

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' SGD TEST\n');
fprintf('=============================================================\n');


%% ============================================================
% TEST 1
%
% One single batch:
%
%       SGD must be identical to full-batch gradient descent.
%
% Since nTrain = 150 < batchSize = 200,
% LossFunctional creates exactly one batch.
%
% Therefore, for every epoch:
%
%   theta_(k+1) = theta_k - lr * grad J(theta_k)
%
% must be exactly the same in both implementations.
% ============================================================

fprintf('\n=== TEST 1: one batch, SGD vs gradient descent ===\n');

rng(0);

lr  = 0.05;
nEp = 20;


%
% No validation set in this first test.
%
d = makeData( ...
    150, ...
    0, ...
    nIn, ...
    nOut);


%
% Small non-zero lambda also checks that the
% regularization term is included consistently.
%
lambda = 1e-3;


[C,lv,~] = ...
    makeCost(d,lambda);


theta0 = ...
    lv.thetavec;


%
% Run SGD.
%
s = struct( ...
    'costFunc',       C, ...
    'designVariable', lv, ...
    'plotter',        [], ...
    'learningRate',   lr, ...
    'maxEpochs',      nEp);


opt = ...
    SGD(s);


opt.compute();


thetaSGD = ...
    lv.thetavec;


%
% Independent reference:
%
% plain full-batch gradient descent.
%
thetaGD = ...
    theta0;


for k = 1:nEp

    C.computeFunctionAndGradient(thetaGD);

    grad = ...
        C.gradient;

    thetaGD = ...
        thetaGD - lr*grad;

end


difference = ...
    norm(thetaSGD(:)-thetaGD(:),2);


relativeDifference = ...
    difference ...
    / ...
    max(norm(thetaGD(:),2),eps);


fprintf('\n');
fprintf( ...
    '||theta_SGD - theta_GD||        = %.16e\n', ...
    difference);

fprintf( ...
    'relative parameter difference   = %.16e\n', ...
    relativeDifference);


if relativeDifference < 1e-12

    fprintf('TEST 1 RESULT: PASS\n');

elseif relativeDifference < 1e-9

    fprintf('TEST 1 RESULT: ACCEPTABLE\n');

else

    fprintf('TEST 1 RESULT: FAIL\n');

end



%% ============================================================
% TEST 2
%
% Early stopping and restoration of the best validation model.
%
% We want to verify:
%
%   1) a finite best validation error is found;
%   2) bestEpoch is recorded;
%   3) after SGD terminates, theta is restored to bestTheta;
%   4) therefore:
%
%       Jval(theta_final) = bestValidationError.
%
% ============================================================

fprintf('\n');
fprintf('=== TEST 2: early stopping / best-theta restoration ===\n');

rng(1);


d = makeData( ...
    600, ...
    200, ...
    nIn, ...
    nOut);


%
% No regularization here so that this test isolates
% validation / early stopping behavior.
%
[C,lv,L] = ...
    makeCost(d,0);


s = struct( ...
    'costFunc',       C, ...
    'designVariable', lv, ...
    'plotter',        [], ...
    'learningRate',   0.2, ...
    'maxEpochs',      300, ...
    'earlyStop',      10);


opt = ...
    SGD(s);


opt.compute();


bestEpoch = ...
    opt.getBestEpoch();


bestValidation = ...
    opt.getBestValidationError();


%
% LossFunctional uses the CURRENT theta stored in lv.
% SGD should already have restored bestTheta.
%
finalValidation = ...
    L.getTestError();


validationDifference = ...
    finalValidation - bestValidation;


fprintf('\n');
fprintf( ...
    'Best epoch                  = %d\n', ...
    bestEpoch);

fprintf( ...
    'Best validation loss        = %.16e\n', ...
    bestValidation);

fprintf( ...
    'Final restored val. loss    = %.16e\n', ...
    finalValidation);

fprintf( ...
    'Final validation - best     = %.16e\n', ...
    validationDifference);


if ...
        bestEpoch > 0 ...
        && ...
        isfinite(bestValidation) ...
        && ...
        abs(validationDifference) < 1e-12

    fprintf('TEST 2 RESULT: PASS\n');

elseif ...
        bestEpoch > 0 ...
        && ...
        isfinite(bestValidation) ...
        && ...
        abs(validationDifference) < 1e-9

    fprintf('TEST 2 RESULT: ACCEPTABLE\n');

else

    fprintf('TEST 2 RESULT: FAIL\n');

end


fprintf('\n');
fprintf('=============================================================\n');
fprintf(' END SGD TEST\n');
fprintf('=============================================================\n');

end



%% ============================================================
% SYNTHETIC DATA
% ============================================================

function d = makeData(nTr,nVal,nIn,nOut)

%
% This test is explicitly written for 2 inputs and 6 outputs.
%
if nIn ~= 2 || nOut ~= 6

    error( ...
        'testSGD:Dimensions', ...
        'This synthetic test expects nIn = 2 and nOut = 6.');

end


f = @(X) [ ...
    sin(3*X(:,1)), ...
    cos(2*X(:,2)), ...
    X(:,1).*X(:,2), ...
    tanh(X(:,1)), ...
    X(:,2).^2, ...
    sin(X(:,1)+X(:,2))];


%% Training data

d.Xtrain = ...
    randn(nTr,nIn);

d.Ytrain = ...
    f(d.Xtrain) ...
    + ...
    0.3*randn(nTr,nOut);


%% Validation data

d.Xtest = ...
    randn(nVal,nIn);

d.Ytest = ...
    f(d.Xtest);


%% Network metadata

d.nFeatures = ...
    nIn;

d.nLabels = ...
    nOut;

end



%% ============================================================
% CREATE NETWORK + LOSS + REGULARIZATION + COST
% ============================================================

function [C,lv,L] = makeCost(d,lambda)

%% Network

sN.hiddenLayers = ...
    [32 32];

sN.HUtype = ...
    'tanh';

sN.OUtype = ...
    'linear';

sN.data = ...
    d;


net = ...
    Network(sN);


%% Learnable variables

lv = ...
    net.getLearnableVariables();


%% Data loss

sL.costType = ...
    'L2';

sL.designVariable = ...
    lv;

sL.network = ...
    net;

sL.data = ...
    d;


L = ...
    LossFunctional(sL);


%% L2 regularization

sR.designVariable = ...
    lv;


R = ...
    Sh_Func_L2norm(sR);


%% Total cost

sC.shapeFunctions = ...
    {L,R};

sC.weights = ...
    [1,lambda];


C = ...
    CostNN(sC);

end