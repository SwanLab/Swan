function s = latticeNNSettings()
    s.seed = 1;
    s.reference.trainFileName = 'latticeTrain.mat';
    s.reference.testFileName  = 'latticeTest.mat';
    s.output.type          = 'cholesky';
    s.output.densityFactor = true;
    s.data.trainFileName   = 'latticeTrainCholesky.mat';
    s.data.testFileName    = 'latticeTestCholesky.mat';
    s.data.xFeatures       = 1:3;
    s.data.yFeatures       = 4:9;
    s.data.polynomialOrder = 1;
    s.networkParams.hiddenLayers = [64 64 64];
    s.networkParams.HUtype       = 'SoftPlus';
    s.networkParams.OUtype       = 'linear';
    s.optimizerParams.type               = 'Adam';
    s.optimizerParams.learningRate       = 1e-3;
    s.optimizerParams.maxEpochs          = 100000;
    s.optimizerParams.earlyStop          = 5000;
    s.optimizerParams.targetCost         = 0;
    s.optimizerParams.restoreBestWeights = true;
    s.optimizerParams.displayEvery       = 1000;
    s.costParams.costType = 'L2';
    s.costParams.lambda   = 1e-9;
    s.parameterNames = {'rho','a','b'};
    s.evaluation.bins = [ ...
        struct('name','rho','map',@(X) X(:,1),'edges',0:0.2:1), ...
        struct('name','|a|','map',@(X) abs(X(:,2)),'edges',0:0.25:0.75), ...
        struct('name','max(b,1/b)','map',@(X) max(X(:,3),1./X(:,3)),'edges',[1 1.5 2 3])];
    s.outputFile     = 'latticeNN_Adam_64x3_cholesky.mat';
end