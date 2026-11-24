function s = latticeSidesNNSettings()
s = latticeNNSettings();
s.reference.trainFileName = 'latticeSidesTrain.mat';
s.reference.testFileName  = 'latticeSidesTest.mat';
s.output.type             = 'cholesky';
s.output.densityFactor    = struct( ...
    'fun',      @(X) 1 - X(:,1).*X(:,2), ...
    'gradient', @(X) [-X(:,2), -X(:,1), zeros(size(X,1),1)]);
s.data.trainFileName      = 'latticeSidesTrainCholesky.mat';
s.data.testFileName       = 'latticeSidesTestCholesky.mat';
s.data.xFeatures          = 1:3;
s.data.yFeatures          = 5:10;
s.transformation          = HoleSidesTransformation();
s.symmetryColumns         = 1:3;
s.parameterNames          = {'m1','m2','a'};
s.evaluation.bins = [ ...
    struct('name','m1', 'map',@(X) X(:,1),     'edges',linspace(0,1,6)), ...
    struct('name','m2', 'map',@(X) X(:,2),     'edges',linspace(0,1,6)), ...
    struct('name','|a|','map',@(X) abs(X(:,3)),'edges',0:0.2:0.8)];
s.outputFile = 'latticeSidesNN_native_Adam_64x3_cholesky.mat';
end