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