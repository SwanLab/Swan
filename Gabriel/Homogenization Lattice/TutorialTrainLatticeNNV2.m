clc
clear
close all

rng(1)

fprintf('\n')
fprintf('============================================\n')
fprintf('          LATTICE NN TRAINING\n')
fprintf('============================================\n')


%% ------------------------------------------------------------
%  DATA
% -------------------------------------------------------------

sData.fileName = ...
    'LatticeNNData.csv';

sData.polynomialOrder = ...
    3;

sData.testRatio = ...
    30;

sData.xFeatures = ...
    1:3;

sData.yFeatures = ...
    4:9;

data = ...
    Data(sData);


fprintf('\nData loaded\n')

fprintf(...
    'Training samples = %d\n',...
    size(data.Xtrain,1))

fprintf(...
    'Test samples     = %d\n',...
    size(data.Xtest,1))

fprintf(...
    'Input features   = %d\n',...
    data.nFeatures)

fprintf(...
    'Output labels    = %d\n',...
    data.nLabels)


%% ------------------------------------------------------------
%  NETWORK PARAMETERS
% -------------------------------------------------------------

s.networkParams.hiddenLayers = ...
    [100 200 300 200 100 50];

s.networkParams.HUtype = ...
    'SoftPlus';

s.networkParams.OUtype = ...
    'linear';


fprintf('\nNetwork architecture\n')

disp(...
    [...
    data.nFeatures,...
    s.networkParams.hiddenLayers,...
    data.nLabels])


%% ------------------------------------------------------------
%  COST PARAMETERS
% -------------------------------------------------------------

s.costParams.costType = ...
    'L2';

s.costParams.lambda = ...
    1e-9;


%% ------------------------------------------------------------
%  OPTIMIZER PARAMETERS
% -------------------------------------------------------------

s.optimizerParams.learningRate = ...
    5e-3;

s.optimizerParams.maxEpochs = ...
    10000;

s.optimizerParams.earlyStop = ...
    1000;


%% ------------------------------------------------------------
%  CREATE OPTIMIZATION PROBLEM
% -------------------------------------------------------------

s.data = ...
    data;

problemNN = ...
    OptimizationProblemNN(s);


fprintf('\n')
fprintf('Starting training...\n')
fprintf('\n')


%% ------------------------------------------------------------
%  TRAIN NETWORK
% -------------------------------------------------------------

problemNN.solve();


%% ------------------------------------------------------------
%  COST HISTORY
% -------------------------------------------------------------

problemNN.plotCostFnc();


%% ------------------------------------------------------------
%  PHYSICAL PREDICTIONS
% -------------------------------------------------------------

YtrainPred = ...
    problemNN.computeOutputValues(...
    data.XtrainPhysical);

YtestPred = ...
    problemNN.computeOutputValues(...
    data.XtestPhysical);


%% ------------------------------------------------------------
%  GLOBAL ERRORS
% -------------------------------------------------------------

trainError = ...
    problemNN.computeError(...
    data.XtrainPhysical,...
    data.YtrainPhysical);

testError = ...
    problemNN.computeError(...
    data.XtestPhysical,...
    data.YtestPhysical);


fprintf('\n')
fprintf('============================================\n')
fprintf('             GLOBAL ERRORS\n')
fprintf('============================================\n')

fprintf(...
    'Training relative error = %.12e\n',...
    trainError)

fprintf(...
    'Test relative error     = %.12e\n',...
    testError)


%% ------------------------------------------------------------
%  ERROR BY CONSTITUTIVE COMPONENT
% -------------------------------------------------------------

componentNames = {...
    'C1111',...
    'C2222',...
    'C1122',...
    'C1112',...
    'C2212',...
    'C1212'};


fprintf('\n')
fprintf('============================================\n')
fprintf('        ERROR BY TENSOR COMPONENT\n')
fprintf('============================================\n')


componentError = ...
    zeros(1,data.nLabels);

componentMAE = ...
    zeros(1,data.nLabels);

componentMaxError = ...
    zeros(1,data.nLabels);


for iC = 1:data.nLabels

    Yref = ...
        data.YtestPhysical(:,iC);

    Ypred = ...
        YtestPred(:,iC);

    componentError(iC) = ...
        norm(Ypred-Yref)/...
        max(norm(Yref),eps);

    componentMAE(iC) = ...
        mean(...
        abs(Ypred-Yref));

    componentMaxError(iC) = ...
        max(...
        abs(Ypred-Yref));

    fprintf(...
        '%s\n',...
        componentNames{iC})

    fprintf(...
        '  Relative error = %.12e\n',...
        componentError(iC))

    fprintf(...
        '  MAE            = %.12e\n',...
        componentMAE(iC))

    fprintf(...
        '  Max abs error  = %.12e\n',...
        componentMaxError(iC))

end


%% ------------------------------------------------------------
%  JACOBIAN CHECK
% -------------------------------------------------------------

nCheck = ...
    min(10,size(data.XtestPhysical,1));

Xcheck = ...
    data.XtestPhysical(1:nCheck,:);

[Ycheck,dYcheck] = ...
    problemNN.computeOutputValuesAndGradient(...
    Xcheck);


fprintf('\n')
fprintf('============================================\n')
fprintf('             JACOBIAN CHECK\n')
fprintf('============================================\n')

fprintf(...
    'Y size  = %d x %d\n',...
    size(Ycheck,1),...
    size(Ycheck,2))

fprintf(...
    'dY size = %d x %d x %d\n',...
    size(dYcheck,1),...
    size(dYcheck,2),...
    size(dYcheck,3))


if any(~isfinite(Ycheck(:)))

    error(...
        'Non-finite NN outputs were detected.')

end


if any(~isfinite(dYcheck(:)))

    error(...
        'Non-finite NN derivatives were detected.')

end


fprintf('\nAll NN outputs are finite.\n')
fprintf('All NN derivatives are finite.\n')


%% ------------------------------------------------------------
%  DERIVATIVE MAGNITUDES
% -------------------------------------------------------------

derivativeNames = {...
    'rho',...
    'a',...
    'b'};


fprintf('\n')
fprintf('============================================\n')
fprintf('        NN DERIVATIVE MAGNITUDES\n')
fprintf('============================================\n')


for iX = 1:data.nFeatures

    derivativeI = ...
        reshape(...
        dYcheck(iX,:,:),...
        data.nLabels,...
        nCheck)';

    fprintf(...
        'dC/d%s: mean Frobenius/sample = %.12e\n',...
        derivativeNames{iX},...
        mean(...
        sqrt(...
        sum(derivativeI.^2,2))))

end


%% ------------------------------------------------------------
%  SAVE MODEL
% -------------------------------------------------------------

save(...
    'LatticeNN_FirstModel.mat',...
    'problemNN',...
    'data',...
    's',...
    'trainError',...
    'testError',...
    'componentError',...
    'componentMAE',...
    'componentMaxError',...
    '-v7.3');


fprintf('\n')
fprintf('Model saved: LatticeNN_FirstModel.mat\n')

fprintf('\n')
fprintf('============================================\n')
fprintf('        LATTICE NN TRAINING COMPLETED\n')
fprintf('============================================\n')
%% ------------------------------------------------------------
%  ERROR BY RHO REGION
% -------------------------------------------------------------

rhoTest = ...
    data.XtestPhysical(:,1);

regions = [...
    0.00 0.20;
    0.20 0.40;
    0.40 0.60;
    0.60 0.80;
    0.80 1.00];


fprintf('\n')
fprintf('============================================\n')
fprintf('             ERROR BY RHO REGION\n')
fprintf('============================================\n')


for iR = 1:size(regions,1)

    rhoMin = ...
        regions(iR,1);

    rhoMax = ...
        regions(iR,2);

    if iR < size(regions,1)

        indices = ...
            rhoTest >= rhoMin & ...
            rhoTest < rhoMax;

    else

        indices = ...
            rhoTest >= rhoMin & ...
            rhoTest <= rhoMax;

    end

    Yref = ...
        data.YtestPhysical(indices,:);

    Ypred = ...
        YtestPred(indices,:);

    errorRho = ...
        norm(...
        Ypred-Yref,...
        'fro')/...
        max(...
        norm(Yref,'fro'),...
        eps);

    fprintf(...
        'rho = [%.2f, %.2f]  -> error = %.12e   N = %d\n',...
        rhoMin,...
        rhoMax,...
        errorRho,...
        sum(indices));

end