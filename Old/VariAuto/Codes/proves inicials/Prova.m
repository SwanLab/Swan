hiddenlayers    = [1,2];

data      = Data('../Datasets/Iris.csv',30,1);
structure = [data.nFeatures,hiddenlayers,data.nLabels];
actual   = optimizationProblem(data,structure);