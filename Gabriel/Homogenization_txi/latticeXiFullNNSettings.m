function s = latticeXiFullNNSettings()
s = latticeNNSettings();
s.reference.trainFileName = 'latticeXiTrain.mat';
s.reference.testFileName  = 'latticeXiTest.mat';
s.data.trainFileName      = 'latticeXiTrainCholesky.mat';
s.data.testFileName       = 'latticeXiTestCholesky.mat';
s.data.xFeatures          = 1:3;
s.data.yFeatures          = 4:9;
s.parameterNames          = {'rho','a','xi'};
s.evaluation.bins = [ ...
    struct('name','rho','map',@(X) X(:,1),'edges',0:0.2:1), ...
    struct('name','|a|','map',@(X) abs(X(:,2)),'edges',0:0.2:0.8), ...
    struct('name','xi','map',@(X) X(:,3),'edges',linspace(0,pi/2,5))];
s.outputFile = 'latticeXiNN_Adam_64x3_cholesky.mat';
end