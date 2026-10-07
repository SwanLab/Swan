function testNetworkInputJacobian()

clc;

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' NETWORK INPUT-DERIVATIVE TEST\n');
fprintf('=============================================================\n');


%% ============================================================
% LOAD DATA
%
% We only use the stored input data.
% The old trained network is NOT used.
% ============================================================

if ~exist('HomogNN.mat','file')

    error( ...
        'testNetworkInputJacobian:MissingFile', ...
        'HomogNN.mat was not found.');

end


S = ...
    load('HomogNN.mat','data');

data = ...
    S.data;


fprintf( ...
    'Training samples : %d\n', ...
    size(data.Xtrain,1));

fprintf( ...
    'Network inputs   : %d\n', ...
    size(data.Xtrain,2));

fprintf( ...
    'Network outputs  : %d\n', ...
    size(data.Ytrain,2));


%% ============================================================
% TEST BOTH ACTIVATIONS
%
% tanh:
%   smooth reference test
%
% ReLU:
%   verifies the activation currently used by the old baseline
% ============================================================

resultsTanh = ...
    runCheck(data,'tanh');

resultsReLU = ...
    runCheck(data,'ReLU');


%% ============================================================
% FINAL SUMMARY
% ============================================================

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' FINAL SUMMARY\n');
fprintf('=============================================================\n');

fprintf( ...
    'tanh : internal Jacobian/directional error = %.3e\n', ...
    resultsTanh.internalError);

fprintf( ...
    'tanh : best FD relative error             = %.3e\n', ...
    resultsTanh.bestFDError);

fprintf( ...
    'tanh : coordinate FD spot-check error     = %.3e\n', ...
    resultsTanh.spotError);

fprintf('\n');

fprintf( ...
    'ReLU : internal Jacobian/directional error = %.3e\n', ...
    resultsReLU.internalError);

fprintf( ...
    'ReLU : best FD relative error             = %.3e\n', ...
    resultsReLU.bestFDError);

fprintf( ...
    'ReLU : coordinate FD spot-check error     = %.3e\n', ...
    resultsReLU.spotError);


%
% The smooth tanh test is the strongest mathematical check.
%
passTanh = ...
    resultsTanh.internalError < 1e-10 ...
    && ...
    resultsTanh.bestFDError < 1e-7 ...
    && ...
    resultsTanh.spotError < 1e-7;


%
% ReLU receives a looser FD tolerance because finite
% differences can cross activation kinks.
%
passReLU = ...
    resultsReLU.internalError < 1e-10 ...
    && ...
    resultsReLU.bestFDError < 1e-4;


fprintf('\n');

if passTanh && passReLU

    fprintf('OVERALL RESULT: PASS\n');

elseif passTanh

    fprintf( ...
        ['OVERALL RESULT: TANH PASS; ', ...
         'inspect ReLU finite-difference behavior.\n']);

else

    fprintf( ...
        ['OVERALL RESULT: FAIL - input derivatives ', ...
         'must be investigated.\n']);

end


fprintf('=============================================================\n');

end



%% ============================================================
% RUN CHECK FOR ONE ACTIVATION
% ============================================================

function results = runCheck(data,activation)

fprintf('\n');
fprintf('-------------------------------------------------------------\n');
fprintf(' ACTIVATION: %s\n',activation);
fprintf('-------------------------------------------------------------\n');


%% ============================================================
% CREATE FRESH NETWORK
%
% Reset the seed so tanh and ReLU start from exactly the same
% initial weights and biases.
% ============================================================

rng(1234);


sN.hiddenLayers = ...
    [150 200 300 200 150 50];

sN.HUtype = ...
    activation;

sN.OUtype = ...
    'linear';

sN.data = ...
    data;


net = ...
    Network(sN);


%% ============================================================
% SELECT A SMALL SET OF INPUT POINTS
%
% We do not need all 2977 points for a derivative check.
% ============================================================

nPtsTest = ...
    min(8,size(data.Xtrain,1));


rng(4321);

idx = ...
    randperm(size(data.Xtrain,1),nPtsTest);


X = ...
    data.Xtrain(idx,:);


[nPts,nIn] = ...
    size(X);


Y0 = ...
    net.computeYOut(X);


nOut = ...
    size(Y0,2);


fprintf( ...
    'Points used      : %d\n', ...
    nPts);

fprintf( ...
    'Input dimension  : %d\n', ...
    nIn);

fprintf( ...
    'Output dimension : %d\n', ...
    nOut);


%% ============================================================
% 1. FULL JACOBIAN
%
% Expected shape:
%
%       nPts x nOut x nIn
%
% J(i,m,q) = dY_m(i) / dX_q(i)
% ============================================================

J = ...
    net.networkJacobian(X);


expectedSize = ...
    [nPts,nOut,nIn];


actualSize = ...
    size(J);


%
% MATLAB can omit a trailing singleton dimension,
% but here nIn = 27, so all three dimensions must exist.
%
if ...
        numel(actualSize) ~= 3 ...
        || ...
        any(actualSize ~= expectedSize)

    error( ...
        'testNetworkInputJacobian:JacobianSize', ...
        ['networkJacobian returned size %s; ', ...
         'expected %s.'], ...
        mat2str(actualSize), ...
        mat2str(expectedSize));

end


if any(~isfinite(J(:)))

    error( ...
        'testNetworkInputJacobian:NonFiniteJacobian', ...
        'networkJacobian contains NaN or Inf.');

end


%% ============================================================
% RANDOM INPUT DIRECTION
%
% Every test point receives its own perturbation direction.
% ============================================================

rng(999);

V = ...
    randn(nPts,nIn);


%
% Normalize the complete perturbation.
%
V = ...
    V / norm(V(:),2);


%% ============================================================
% 2. CONTRACT FULL JACOBIAN WITH V
%
%       dY = J_X V
%
% Result:
%
%       nPts x nOut
% ============================================================

dYfromJacobian = ...
    zeros(nPts,nOut);


for q = 1:nIn

    dYfromJacobian = ...
        dYfromJacobian ...
        + ...
        J(:,:,q).*V(:,q);

end


%% ============================================================
% 3. NETWORK DIRECTIONAL DERIVATIVE
% ============================================================

dX = ...
    zeros(nPts,nIn,1);

dX(:,:,1) = ...
    V;


dYdirectional = ...
    net.networkDirectionalDerivative(X,dX);


%
% Accept both
%
%   nPts x nOut
%
% and
%
%   nPts x nOut x 1.
%
if ndims(dYdirectional) == 3

    dYdirectional = ...
        dYdirectional(:,:,1);

end


if ~isequal(size(dYdirectional),[nPts,nOut])

    error( ...
        'testNetworkInputJacobian:DirectionalSize', ...
        ['networkDirectionalDerivative returned size %s; ', ...
         'expected [%d %d].'], ...
        mat2str(size(dYdirectional)), ...
        nPts,nOut);

end


%% ============================================================
% INTERNAL CONSISTENCY
%
% networkJacobian and networkDirectionalDerivative should
% represent the same derivative.
% ============================================================

internalAbs = ...
    norm( ...
        dYfromJacobian(:) ...
        - ...
        dYdirectional(:), ...
        2);


internalDen = ...
    max( ...
        norm(dYfromJacobian(:),2), ...
        eps);


internalError = ...
    internalAbs/internalDen;


fprintf('\n');
fprintf('1) Jacobian vs directional derivative\n');

fprintf( ...
    '   relative error = %.16e\n', ...
    internalError);


%% ============================================================
% 4. DIRECTIONAL DERIVATIVE VS CENTRAL FINITE DIFFERENCE
%
%       Y(X+hV)-Y(X-hV)
%       ----------------
%              2h
%
% ============================================================

hValues = ...
    [1e-2 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8];


fdErrors = ...
    zeros(size(hValues));


fprintf('\n');
fprintf('2) Directional derivative vs finite differences\n');
fprintf('\n');

fprintf( ...
    ['----------------------------------------------------------------', ...
     '-------------\n']);

fprintf( ...
    '       h              abs error             relative error\n');

fprintf( ...
    ['----------------------------------------------------------------', ...
     '-------------\n']);


for ih = 1:numel(hValues)

    h = ...
        hValues(ih);


    Xplus = ...
        X + h*V;

    Xminus = ...
        X - h*V;


    Yplus = ...
        net.computeYOut(Xplus);

    Yminus = ...
        net.computeYOut(Xminus);


    dYfd = ...
        (Yplus-Yminus)/(2*h);


    absError = ...
        norm( ...
            dYfd(:) ...
            - ...
            dYfromJacobian(:), ...
            2);


    denom = ...
        max( ...
            norm(dYfromJacobian(:),2), ...
            eps);


    relError = ...
        absError/denom;


    fdErrors(ih) = ...
        relError;


    fprintf( ...
        '%12.1e    %20.12e    %20.12e\n', ...
        h, ...
        absError, ...
        relError);

end


[bestFDError,idBest] = ...
    min(fdErrors);


fprintf( ...
    ['----------------------------------------------------------------', ...
     '-------------\n']);

fprintf( ...
    'Best FD relative error = %.6e at h = %.1e\n', ...
    bestFDError, ...
    hValues(idBest));


%% ============================================================
% 5. SPOT CHECK INDIVIDUAL JACOBIAN ENTRIES
%
% Perturb one input coordinate of one point at a time:
%
%   dY_m / dX_q
%
% and compare the complete output vector against J(i,:,q).
%
% This is independent of the random-direction contraction.
% ============================================================

fprintf('\n');
fprintf('3) Individual Jacobian coordinate spot checks\n');


%
% h = 1e-5 is normally a good compromise for tanh.
% For ReLU, a sampled point may occasionally cross a kink.
%
hSpot = ...
    1e-5;


nChecks = ...
    min(12,nPts*nIn);


rng(2026);


pointIndex = ...
    randi(nPts,[nChecks,1]);

inputIndex = ...
    randi(nIn,[nChecks,1]);


spotNumerator = ...
    0;

spotDenominator = ...
    0;


for k = 1:nChecks

    ip = ...
        pointIndex(k);

    iq = ...
        inputIndex(k);


    Xplus = ...
        X;

    Xminus = ...
        X;


    Xplus(ip,iq) = ...
        Xplus(ip,iq) + hSpot;

    Xminus(ip,iq) = ...
        Xminus(ip,iq) - hSpot;


    Yplus = ...
        net.computeYOut(Xplus);

    Yminus = ...
        net.computeYOut(Xminus);


    %
    % Only output ip depends on input ip because network
    % samples are evaluated independently.
    %
    dYfd = ...
        (Yplus(ip,:)-Yminus(ip,:))/(2*hSpot);


    dYexact = ...
        reshape(J(ip,:,iq),1,nOut);


    spotNumerator = ...
        spotNumerator ...
        + ...
        sum((dYfd-dYexact).^2);


    spotDenominator = ...
        spotDenominator ...
        + ...
        sum(dYexact.^2);

end


spotError = ...
    sqrt(spotNumerator) ...
    / ...
    max(sqrt(spotDenominator),eps);


fprintf( ...
    '   checks performed = %d\n', ...
    nChecks);

fprintf( ...
    '   h                 = %.1e\n', ...
    hSpot);

fprintf( ...
    '   relative error    = %.16e\n', ...
    spotError);


%% ============================================================
% RESULTS
% ============================================================

results.internalError = ...
    internalError;

results.bestFDError = ...
    bestFDError;

results.bestH = ...
    hValues(idBest);

results.spotError = ...
    spotError;


fprintf('\n');


if strcmpi(activation,'tanh')

    pass = ...
        internalError < 1e-10 ...
        && ...
        bestFDError < 1e-7 ...
        && ...
        spotError < 1e-7;

else

    %
    % ReLU finite differences may cross activation boundaries.
    %
    pass = ...
        internalError < 1e-10 ...
        && ...
        bestFDError < 1e-4;

end


if pass

    fprintf('ACTIVATION RESULT: PASS\n');

else

    fprintf('ACTIVATION RESULT: INSPECT\n');

end

end