clc
clear
close all

rng(1)

lattice = TutorialLatticeNN();

rho = lattice.rho(:);
a = lattice.a(:);
b = lattice.b(:);
C = lattice.Chomog;

nRho = length(rho);
nA = length(a);
nB = length(b);

Cminor1 = permute(C,[2 1 3 4 5 6 7]);
Cminor2 = permute(C,[1 2 4 3 5 6 7]);
Cmajor = permute(C,[3 4 1 2 5 6 7]);

errMinor1 = max(abs(C(:)-Cminor1(:)));
errMinor2 = max(abs(C(:)-Cminor2(:)));
errMajor = max(abs(C(:)-Cmajor(:)));

fprintf('Minor symmetry ij: %.12e\n',errMinor1)
fprintf('Minor symmetry kl: %.12e\n',errMinor2)
fprintf('Major symmetry:    %.12e\n',errMajor)

[R,A,B] = ndgrid(rho,a,b);

C1111 = reshape(C(1,1,1,1,:,:,:),nRho,nA,nB);
C2222 = reshape(C(2,2,2,2,:,:,:),nRho,nA,nB);
C1122 = reshape(C(1,1,2,2,:,:,:),nRho,nA,nB);
C1112 = reshape(C(1,1,1,2,:,:,:),nRho,nA,nB);
C2212 = reshape(C(2,2,1,2,:,:,:),nRho,nA,nB);
C1212 = reshape(C(1,2,1,2,:,:,:),nRho,nA,nB);

dataset = [R(:),A(:),B(:),...
           C1111(:),C2222(:),C1122(:),...
           C1112(:),C2212(:),C1212(:)];

X = dataset(:,1:3);
Y = dataset(:,4:9);

disp(' ')
disp('Physical input statistics')
disp(' ')
disp('              min           max           mean          std')
disp([min(X,[],1)', ...
      max(X,[],1)', ...
      mean(X,1)', ...
      std(X,0,1)'])

disp(' ')
disp('Physical output statistics')
disp(' ')
disp('              min           max           mean          std')
disp([min(Y,[],1)', ...
      max(Y,[],1)', ...
      mean(Y,1)', ...
      std(Y,0,1)'])

disp(' ')
disp('Maximum absolute value per physical output')
disp(max(abs(Y),[],1))

fileName = 'LatticeNNData.csv';

writematrix(dataset,fileName)

fprintf('\nNumber of samples: %d\n',size(dataset,1))
fprintf('Number of inputs:  %d\n',3)
fprintf('Number of outputs: %d\n\n',6)

s.fileName = fileName;

s.polynomialOrder = 3;
s.testRatio = 30;

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
disp('Maximum absolute normalized training output')
disp(max(abs(data.Ytrain),[],1))

s.data = data;

T_NN = OptimizationProblemNN(s);

T_NN.solve();

T_NN.plotCostFnc();

Ytrain = T_NN.computeOutputValues(...
    data.XtrainPhysical);

Ytest = T_NN.computeOutputValues(...
    data.XtestPhysical);

errorTrain = T_NN.computeError(...
    data.XtrainPhysical,...
    data.YtrainPhysical);

errorTest = T_NN.computeError(...
    data.XtestPhysical,...
    data.YtestPhysical);

fprintf('\nRelative training error = %.12e\n',errorTrain)
fprintf('Relative test error     = %.12e\n',errorTest)

errorComponentsTrain = ...
    vecnorm(Ytrain-data.YtrainPhysical,2,1)./...
    max(vecnorm(data.YtrainPhysical,2,1),eps);

errorComponentsTest = ...
    vecnorm(Ytest-data.YtestPhysical,2,1)./...
    max(vecnorm(data.YtestPhysical,2,1),eps);

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