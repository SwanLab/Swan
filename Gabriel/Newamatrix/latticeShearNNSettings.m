function s = latticeShearNNSettings()
s = latticeNNSettings();
s.reference.trainFileName = 'latticeShearTrain.mat';
s.reference.testFileName  = 'latticeShearTest.mat';
s.data.trainFileName      = 'latticeShearTrainCholesky.mat';
s.data.testFileName       = 'latticeShearTestCholesky.mat';
s.data.xFeatures          = 1:3;
s.data.yFeatures          = 4:9;
s.parameterNames          = {'rho','a','xi'};
s.transformation          = ShearTransformation();
s.evaluation.bins = [ ...
    struct('name','rho','map',@(X) X(:,1),'edges',0:0.2:1), ...
    struct('name','|a|','map',@(X) abs(X(:,2)),'edges',0:0.5:2), ...
    struct('name','xi','map',@(X) X(:,3),'edges',linspace(0,pi/2,5))];
s.outputFile = 'latticeShearNN_Adam_64x3_cholesky.mat';
end