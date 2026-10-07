clear;
clc;

rng(0);

%% ============================================================
%  LOAD TRAINING DATA
% ============================================================

if ~exist('HomogNN.mat','file')
    error(['HomogNN.mat was not found. ', ...
           'Run the data preparation/fitting once before this test.']);
end

S = load('HomogNN.mat','data');

data = S.data;


%% ============================================================
%  CREATE A FRESH NETWORK
%
%  IMPORTANT:
%  We do NOT train the network here.
%  The purpose is only to verify:
%
%       analytical dJ/dtheta
%
%  against
%
%       finite differences of J(theta).
%
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


fprintf('\n=============================================\n');
fprintf(' NN TRAINING GRADIENT CHECK\n');
fprintf('=============================================\n');

fprintf('Number of parameters: %d\n',numel(theta0));
fprintf('Number of training samples: %d\n',size(data.Xtrain,1));
fprintf('Number of network inputs: %d\n',size(data.Xtrain,2));
fprintf('Number of outputs: %d\n',size(data.Ytrain,2));


%% ============================================================
%  DATA LOSS
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
%
%  Use lambda = 0 first because this is the current setting
%  of DamageHomogenizationFitter.
% ============================================================

sReg.designVariable = ...
    thetaVar;

regularization = ...
    Sh_Func_L2norm(sReg);

lambda = 0;


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
%  ANALYTICAL GRADIENT
% ============================================================

cost.computeFunctionAndGradient(theta0);

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


fprintf('\nJ(theta) = %.16e\n',J0);
fprintf('||grad J||_2 = %.16e\n',norm(g(:),2));


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
%       dJ/ds = grad(J)^T v
%
dJ_analytical = ...
    sum(g(:).*v(:));


fprintf('\nAnalytical directional derivative:\n');
fprintf('grad(J)^T v = %.16e\n\n',dJ_analytical);


%% ============================================================
%  CENTRAL FINITE DIFFERENCE
%
%          J(theta+h*v) - J(theta-h*v)
% dJ/ds ~= ---------------------------------
%                        2h
%
% ============================================================

hValues = ...
    [1e-2 1e-3 1e-4 1e-5 1e-6 1e-7];

fprintf('-------------------------------------------------------------\n');
fprintf('       h              FD derivative        relative error\n');
fprintf('-------------------------------------------------------------\n');


errors = zeros(size(hValues));


for ih = 1:numel(hValues)

    h = ...
        hValues(ih);


    %% theta + h*v

    thetaPlus = ...
        theta0 + h*v;

    cost.computeFunctionAndGradient(thetaPlus);

    Jplus = ...
        cost.value;


    %% theta - h*v

    thetaMinus = ...
        theta0 - h*v;

    cost.computeFunctionAndGradient(thetaMinus);

    Jminus = ...
        cost.value;


    %% Central finite difference

    dJ_FD = ...
        (Jplus-Jminus)/(2*h);


    %% Relative error

    denom = ...
        max( ...
            [1, ...
             abs(dJ_analytical), ...
             abs(dJ_FD)]);

    errors(ih) = ...
        abs(dJ_FD-dJ_analytical)/denom;


    fprintf( ...
        '%12.1e    %20.12e    %12.4e\n', ...
        h, ...
        dJ_FD, ...
        errors(ih));

end


%% ============================================================
%  RESTORE ORIGINAL PARAMETERS
% ============================================================

thetaVar.thetavec = ...
    theta0;


%% ============================================================
%  RESULT
% ============================================================

[minError,idBest] = ...
    min(errors);

fprintf('-------------------------------------------------------------\n');

fprintf( ...
    '\nBest relative error = %.6e at h = %.1e\n', ...
    minError, ...
    hValues(idBest));


if minError < 1e-6

    fprintf('\nRESULT: PASS\n');

elseif minError < 1e-4

    fprintf('\nRESULT: ACCEPTABLE, but inspect convergence with h.\n');

else

    fprintf('\nRESULT: FAIL - training gradient must be investigated.\n');

end

fprintf('=============================================\n');

function testLossFunctional()

rng(0)
nTr = 500;  nIn = 2;  nOut = 6;  h = 1e-6;
d.nFeatures = nIn;  d.nLabels = nOut;
d.Xtrain = randn(nTr,nIn);   d.Ytrain = randn(nTr,nOut);
d.Xtest  = randn(100,nIn);   d.Ytest  = randn(100,nOut);

sN.hiddenLayers = [8 8];  sN.HUtype = 'tanh';  sN.OUtype = 'linear';  sN.data = d;
net = Network(sN);  lv = net.getLearnableVariables();  th0 = lv.thetavec;

sL.costType = 'L2';  sL.designVariable = lv;  sL.network = net;  sL.data = d;
L = LossFunctional(sL);

fprintf('\n=== LossFunctional ===\n');

% ---- 1) gradiente vs diferencas finitas, mesmo batch (moveBatch = false)
[~,g] = L.computeStochasticCostAndGradient(th0,false);
idx = randperm(numel(th0),30);  gFD = zeros(1,30);
for q = 1:30
    tp = th0;  tp(idx(q)) = tp(idx(q)) + h;
    tm = th0;  tm(idx(q)) = tm(idx(q)) - h;
    jp = L.computeStochasticCostAndGradient(tp,false);
    jm = L.computeStochasticCostAndGradient(tm,false);
    gFD(q) = (jp-jm)/(2*h);
end
fprintf('erro rel. gradiente   = %.2e   (esperado ~1e-9)\n', ...
    norm(g(idx)-gFD)/norm(gFD));

% ---- 2) chamadas por epoca (500 pontos, batch 200 -> 2 batches)
calls = zeros(1,4);  j1 = zeros(1,4);
for ep = 1:4
    [j1(ep),~,isBD] = L.computeStochasticCostAndGradient(th0,true);
    calls(ep) = 1;
    while ~isBD
        [~,~,isBD] = L.computeStochasticCostAndGradient(th0,true);
        calls(ep) = calls(ep) + 1;
    end
end
fprintf('chamadas por epoca    = %s   (esperado [2 2 2 2])\n', mat2str(calls));

% ---- 3) reembaralhamento: 1o batch muda de epoca para epoca
fprintf('custo do 1o batch     = %s\n', mat2str(j1,6));
fprintf('  todos distintos?    = %d   (esperado 1)\n', numel(unique(j1)) == 4);

% ---- 4) getTestError = mesma perda, no conjunto de teste
lv.thetavec = th0;
e = net.computeYOut(d.Xtest) - d.Ytest;
fprintf('getTestError - ref    = %.2e   (esperado 0)\n', ...
    L.getTestError() - 0.5*mean(sum(e.^2,2)));

end

function testSGD()

nIn = 2;  nOut = 6;
fprintf('\n=== SGD ===\n');

% ---- 1) batch unico: SGD == descida do gradiente
rng(0)
lr = 0.05;  nEp = 20;
d = makeData(150,0,nIn,nOut);                  % 150 < 200 -> 1 batch
[C,lv,~] = makeCost(d,1e-3);
th0 = lv.thetavec;

s = struct('costFunc',C,'designVariable',lv,'plotter',[], ...
    'learningRate',lr,'maxEpochs',nEp);
opt = SGD(s);  opt.compute();
thSGD = lv.thetavec;

th = th0;
for k = 1:nEp
    C.computeFunctionAndGradient(th);
    th = th - lr*C.gradient;
end
fprintf('||theta_SGD - theta_GD|| = %.2e   (esperado ~1e-15)\n', norm(thSGD-th));

% ---- 2) early stopping restaura o melhor theta
rng(1)
d = makeData(600,200,nIn,nOut);
[C,lv,L] = makeCost(d,0);
s = struct('costFunc',C,'designVariable',lv,'plotter',[], ...
    'learningRate',0.2,'maxEpochs',300,'earlyStop',10);
opt = SGD(s);  opt.compute();

fprintf('melhor epoca             = %d\n', opt.getBestEpoch());
fprintf('melhor validacao         = %.6e\n', opt.getBestValidationError());
fprintf('testError(final) - melhor = %.2e   (esperado 0)\n', ...
    L.getTestError() - opt.getBestValidationError());

end

function d = makeData(nTr,nTe,nIn,nOut)
f = @(X) [sin(3*X(:,1)), cos(2*X(:,2)), X(:,1).*X(:,2), ...
    tanh(X(:,1)), X(:,2).^2, sin(X(:,1)+X(:,2))];
d.nFeatures = nIn;  d.nLabels = nOut;
d.Xtrain = randn(nTr,nIn);  d.Ytrain = f(d.Xtrain) + 0.3*randn(nTr,nOut);
d.Xtest  = randn(nTe,nIn);  d.Ytest  = f(d.Xtest);
end

function [C,lv,L] = makeCost(d,lambda)
sN = struct('hiddenLayers',[32 32],'HUtype','tanh','OUtype','linear','data',d);
net = Network(sN);  lv = net.getLearnableVariables();
L = LossFunctional(struct('costType','L2','designVariable',lv, ...
    'network',net,'data',d));
R = Sh_Func_L2norm(struct('designVariable',lv));
C = CostNN(struct('shapeFunctions',{{L,R}},'weights',[1 lambda]));
end

% ... mantenha tudo ate a criacao de 'cost', com HUtype = 'tanh' e lambda = 1e-4 ...

cost.setBatchMover(false);          % batch fixo durante todo o teste

cost.computeStochasticFunctionAndGradient(theta0);
J0 = cost.value;
g  = cost.gradient;

if ~isequal(size(g),size(theta0))
    error('Gradient shape mismatch: theta is %s, gradient is %s.', ...
        mat2str(size(theta0)), mat2str(size(g)));
end

v = randn(size(theta0));  v = v/norm(v(:),2);
dJ_analytical = sum(g(:).*v(:));

fprintf('\nJ_batch(theta) = %.16e\n',J0);
fprintf('grad(J)^T v    = %.16e\n\n',dJ_analytical);

hValues = [1e-2 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8];
for ih = 1:numel(hValues)
    h = hValues(ih);
    cost.computeStochasticFunctionAndGradient(theta0 + h*v);  Jp = cost.value;
    cost.computeStochasticFunctionAndGradient(theta0 - h*v);  Jm = cost.value;
    dJ_FD = (Jp-Jm)/(2*h);
    err = abs(dJ_FD-dJ_analytical)/abs(dJ_analytical);     % denominador real
    fprintf('%12.1e    %20.12e    %12.4e\n',h,dJ_FD,err);
end

thetaVar.thetavec = theta0;