function s = latticeSquareTurnNNSettings()
s = latticeNNSettings();
s.reference.trainFileName = 'latticeSquareTurnTrain.mat';
s.reference.testFileName  = 'latticeSquareTurnTest.mat';
s.data.trainFileName      = 'latticeSquareTurnTrainCholesky.mat';
s.data.testFileName       = 'latticeSquareTurnTestCholesky.mat';
s.data.xFeatures          = 1:3;
s.data.yFeatures          = 4:9;
s.parameterNames          = {'rho','a','xi'};
s.transformation          = SquareTurnTransformation();
s.optimizerParams.maxEpochs = 100000;
s.optimizerParams.earlyStop = 5000;
s.evaluation.bins = [ ...
    struct('name','rho','map',@(X) X(:,1),'edges',0:0.2:1), ...
    struct('name','|a|','map',@(X) abs(X(:,2)),'edges',0:0.25:1), ...
    struct('name','xi','map',@(X) X(:,3),'edges',linspace(0,pi/2,5))];
s.outputFile = 'latticeSquareTurnNN_Adam_64x3_cholesky.mat';
end