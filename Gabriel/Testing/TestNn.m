clc
clear
close all

rng(1)

trainFileName = 'LatticeNNTrainData.csv';
testFileName = 'LatticeNNTestData.csv';

trainDataset = readmatrix(trainFileName);
testDataset = readmatrix(testFileName);

XtrainPhysical = trainDataset(:,1:3);
YtrainPhysical = trainDataset(:,4:9);

XtestPhysical = testDataset(:,1:3);
YtestPhysical = testDataset(:,4:9);

disp(' ')
disp('Training physical input statistics')
disp(' ')
disp('              min           max           mean          std')
disp([min(XtrainPhysical,[],1)', ...
      max(XtrainPhysical,[],1)', ...
      mean(XtrainPhysical,1)', ...
      std(XtrainPhysical,0,1)'])

disp(' ')
disp('Test physical input statistics')
disp(' ')
disp('              min           max           mean          std')
disp([min(XtestPhysical,[],1)', ...
      max(XtestPhysical,[],1)', ...
      mean(XtestPhysical,1)', ...
      std(XtestPhysical,0,1)'])

disp(' ')
disp('Training physical output statistics')
disp(' ')
disp('              min           max           mean          std')
disp([min(YtrainPhysical,[],1)', ...
      max(YtrainPhysical,[],1)', ...
      mean(YtrainPhysical,1)', ...
      std(YtrainPhysical,0,1)'])

disp(' ')
disp('Test physical output statistics')
disp(' ')
disp('              min           max           mean          std')
disp([min(YtestPhysical,[],1)', ...
      max(YtestPhysical,[],1)', ...
      mean(YtestPhysical,1)', ...
      std(YtestPhysical,0,1)'])

disp(' ')
disp('Maximum absolute value per training output')
disp(max(abs(YtrainPhysical),[],1))

disp(' ')
disp('Maximum absolute value per test output')
disp(max(abs(YtestPhysical),[],1))

fprintf('\nTraining samples: %d\n',size(trainDataset,1))
fprintf('Test samples:     %d\n',size(testDataset,1))
fprintf('Total samples:    %d\n',...
    size(trainDataset,1)+size(testDataset,1))
fprintf('Number of inputs:  %d\n',3)
fprintf('Number of outputs: %d\n\n',6)

s.trainFileName = trainFileName;
s.testFileName = testFileName;

s.polynomialOrder = 3;

s.networkParams.hiddenLayers = [40 72 150 72];

s.optimizerParams.learningRate = 0.32;
s.optimizerParams.maxEpochs = 1000;

s.costParams.lambda = 0;
s.costParams.costType = 'L2';

s.networkParams.HUtype = 'ReLU';
s.networkParams.OUtype = 'linear';

s.xFeatures = 1:3;
s.yFeatures = 4:9;

data = Data(s);

disp(' ')
disp('Normalized training input statistics')
disp(' ')
disp('              min           max           mean          std')
disp([min(data.Xtrain,[],1)', ...
      max(data.Xtrain,[],1)', ...
      mean(data.Xtrain,1)', ...
      std(data.Xtrain,0,1)'])

disp(' ')
disp('Normalized test input statistics')
disp(' ')
disp('              min           max           mean          std')
disp([min(data.Xtest,[],1)', ...
      max(data.Xtest,[],1)', ...
      mean(data.Xtest,1)', ...
      std(data.Xtest,0,1)'])

disp(' ')
disp('Maximum absolute normalized training output')
disp(max(abs(data.Ytrain),[],1))

disp(' ')
disp('Maximum absolute normalized test output')
disp(max(abs(data.Ytest),[],1))

s.data = data;

T_NN = OptimizationProblemNN(s);

T_NN.solve();

T_NN.plotCostFnc();

Ytrain = ...
    T_NN.computeOutputValues(...
    data.XtrainPhysical);

Ytest = ...
    T_NN.computeOutputValues(...
    data.XtestPhysical);

errorTrain = ...
    T_NN.computeError(...
    data.XtrainPhysical,...
    data.YtrainPhysical);

errorTest = ...
    T_NN.computeError(...
    data.XtestPhysical,...
    data.YtestPhysical);

fprintf('\nRelative training error = %.12e\n',errorTrain)
fprintf('Relative test error     = %.12e\n',errorTest)

errorComponentsTrain = ...
    vecnorm(...
    Ytrain-data.YtrainPhysical,...
    2,1)./...
    max(...
    vecnorm(data.YtrainPhysical,2,1),...
    eps);

errorComponentsTest = ...
    vecnorm(...
    Ytest-data.YtestPhysical,...
    2,1)./...
    max(...
    vecnorm(data.YtestPhysical,2,1),...
    eps);

disp(' ')
disp('Relative training error per component')
disp(' ')
disp(' C1111       C2222       C1122       C1112       C2212       C1212')
disp(errorComponentsTrain)

disp(' ')
disp('Relative test error per component')
disp(' ')
disp(' C1111       C2222       C1122       C1112       C2212       C1212')
disp(errorComponentsTest)

[Ypred,J] = ...
    T_NN.computeOutputValuesAndGradient(...
    data.XtestPhysical);

disp(' ')
disp('Size physical Y')
disp(size(Ypred))

disp('Size physical Jacobian')
disp(size(J))

save('LatticeNN.mat','T_NN')