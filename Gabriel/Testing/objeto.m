clc
clear
close all

rng(1)

trainFileName = ...
    'LatticeNNTrainData.csv';

testFileName = ...
    'LatticeNNTestData.csv';

s.trainFileName = ...
    trainFileName;

s.testFileName = ...
    testFileName;

s.polynomialOrder = ...
    3;

s.networkParams.hiddenLayers = ...
    [40 72 150 72];

s.optimizerParams.learningRate = ...
    0.0005;

s.optimizerParams.maxEpochs = ...
    20000;

s.costParams.lambda = ...
    0;

s.costParams.costType = ...
    'L2';

s.networkParams.HUtype = ...
    'SoftPlus';

s.networkParams.OUtype = ...
    'linear';

s.xFeatures = ...
    1:3;

s.yFeatures = ...
    4:9;

data = ...
    Data(s);

s.data = ...
    data;

aMin = ...
    0;

aMax = ...
    0.98;

bMin = ...
    0.3;

bMax = ...
    3;

hRelative = ...
    5e-3;

s.sensitivityParams.hA = ...
    hRelative*(aMax-aMin);

s.sensitivityParams.hB = ...
    hRelative*(bMax-bMin);

s.sensitivityParams.sensitivityFileName = ...
    'LatticeNNSensitivityData.csv';

s.sensitivityParams.weight = ...
    1;

s.balancingParams.enabled = ...
    true;

s.balancingParams.type = ...
    'ReLoBRaLo';

s.balancingParams.alpha = ...
    0.999;

s.balancingParams.temperature = ...
    0.1;

s.balancingParams.rhoExpected = ...
    0.9999;

s.balancingParams.balanceMask = ...
    [true true false];

T_NN = ...
    OptimizationProblemNN(s);
T_NN.solve();
T_NN.printLossDiagnostics();

h = ...
    T_NN.getBalancingHistory();

fprintf('\n')
fprintf('============================================\n')
fprintf('       RELOBRALO TRAINING SUMMARY\n')
fprintf('============================================\n')

fprintf('\nNumber of balancing updates = %d\n',...
    size(h.weights,2))

fprintf('\nInitial losses\n')
disp(h.loss(:,1)')

fprintf('Final losses\n')
disp(h.loss(:,end)')

fprintf('\nInitial weights\n')
disp(h.weights(:,1)')

fprintf('Final weights\n')
disp(h.weights(:,end)')

fprintf('\nMinimum weights\n')
disp(min(h.weights,[],2)')

fprintf('Maximum weights\n')
disp(max(h.weights,[],2)')

fprintf('Mean weights\n')
disp(mean(h.weights,2)')
figure

plot(...
    1:size(h.loss,2),...
    h.loss(1,:),...
    'LineWidth',1.5)

hold on

plot(...
    1:size(h.loss,2),...
    h.loss(2,:),...
    'LineWidth',1.5)

xlabel('Optimization step')
ylabel('Loss')
legend('L_1','L_2')
grid on


figure

semilogy(...
    1:size(h.loss,2),...
    h.loss(1,:),...
    'LineWidth',1.5)

hold on

semilogy(...
    1:size(h.loss,2),...
    h.loss(2,:),...
    'LineWidth',1.5)

xlabel('Optimization step')
ylabel('Loss')
legend('L_1','L_2')
grid on