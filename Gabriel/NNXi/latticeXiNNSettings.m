function s = latticeXiNNSettings(tag)
if nargin < 1 || isempty(tag)
    tag = 'Rho08';
end
s = latticeNNSettings();
s.reference.trainFileName = ['latticeXi',tag,'Train.mat'];
s.reference.testFileName  = ['latticeXi',tag,'Test.mat'];
s.output.densityFactor    = false;
s.data.trainFileName      = ['latticeXi',tag,'TrainCholesky.mat'];
s.data.testFileName       = ['latticeXi',tag,'TestCholesky.mat'];
s.data.xFeatures          = 2:3;
s.data.yFeatures          = 4:9;
s.parameterNames          = {'a','xi'};
s.evaluation.bins = [ ...
    struct('name','|a|','map',@(X) abs(X(:,1)),'edges',0:0.2:0.8), ...
    struct('name','xi','map',@(X) X(:,2),'edges',linspace(0,pi/2,5))];
s.optimizerParams.maxEpochs = 100000;
s.optimizerParams.earlyStop = 5000;
s.outputFile = ['latticeXiNN_',lower(tag),'_Adam_64x3_cholesky.mat'];
end