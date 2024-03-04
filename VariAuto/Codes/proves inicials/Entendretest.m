%% Initialization of hyperparameters
hiddenlayers    = [1,2];

data      = Data('../Datasets/Iris.csv',30,1);
structure = [data.nFeatures,hiddenlayers,data.nLabels];
network   = Network(data,structure);