clc
clear
close all

rng(1)

s.trainFileName = ...
    'LatticeNNTrainData.csv';

s.testFileName = ...
    'LatticeNNTestData.csv';

s.polynomialOrder = 3;

s.xFeatures = 1:3;
s.yFeatures = 4:9;

data = Data(s);

networkParams.hiddenLayers = ...
    [40 72 150 72];

networkParams.HUtype = ...
    'SoftPlus';

networkParams.OUtype = ...
    'linear';

networkParams.data = data;

network = Network(networkParams);

X = [...
    -0.70 -0.40 -0.20;
    -0.30  0.20  0.60;
     0.10 -0.60  0.30;
     0.50  0.40 -0.70;
     0.80  0.10  0.50];

[Y,J] = ...
    network.computeYOutAndJacobian(X);

h = 1e-6;

nSamples = size(X,1);
nFeatures = size(X,2);
nOutputs = size(Y,2);

Jfd = zeros(...
    nFeatures,...
    nOutputs,...
    nSamples);

for iFeature = 1:nFeatures

    Xp = X;
    Xm = X;

    Xp(:,iFeature) = ...
        Xp(:,iFeature) + h;

    Xm(:,iFeature) = ...
        Xm(:,iFeature) - h;

    Yp = ...
        network.computeYOut(Xp);

    Ym = ...
        network.computeYOut(Xm);

    dY = ...
        (Yp-Ym)/(2*h);

    Jfd(iFeature,:,:) = ...
        permute(dY,[3 2 1]);

end

absError = ...
    max(abs(J(:)-Jfd(:)));

relError = ...
    norm(J(:)-Jfd(:))/...
    max(norm(Jfd(:)),eps);

fprintf('\n')
fprintf('============================================\n')
fprintf('       SOFTPLUS JACOBIAN VALIDATION\n')
fprintf('============================================\n')

fprintf('\nNumber of samples:  %d\n',nSamples)
fprintf('Number of inputs:   %d\n',nFeatures)
fprintf('Number of outputs:  %d\n',nOutputs)

fprintf('\nSize Y\n')
disp(size(Y))

fprintf('Size analytical Jacobian\n')
disp(size(J))

fprintf('Size FD Jacobian\n')
disp(size(Jfd))

fprintf('\nGlobal errors\n')
fprintf('Max absolute error = %.12e\n',absError)
fprintf('Relative error     = %.12e\n',relError)

featureNames = {'rho','a','b'};

fprintf('\n')
fprintf('Error per input feature\n')

for iFeature = 1:nFeatures

    Ji = reshape(...
        J(iFeature,:,:),...
        nOutputs,nSamples);

    Jfdi = reshape(...
        Jfd(iFeature,:,:),...
        nOutputs,nSamples);

    absErrorFeature = ...
        max(abs(Ji(:)-Jfdi(:)));

    relErrorFeature = ...
        norm(Ji(:)-Jfdi(:))/...
        max(norm(Jfdi(:)),eps);

    fprintf(...
        '%s: max abs = %.12e   rel = %.12e\n',...
        featureNames{iFeature},...
        absErrorFeature,...
        relErrorFeature)

end

fprintf('\n')
fprintf('Analytical Jacobian - sample 1\n')
disp(J(:,:,1))

fprintf('Finite difference Jacobian - sample 1\n')
disp(Jfd(:,:,1))

fprintf('Difference - sample 1\n')
disp(J(:,:,1)-Jfd(:,:,1))